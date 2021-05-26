from __future__ import annotations
from typing import Union
from rec import TimebinInterval
from rec_db import *
from rec_format import *
from versioning import *
from dPCA import dPCA
from copy import copy
from itertools import chain
from sklearn.preprocessing import StandardScaler


class Analysis:

    def __init__(self, db: DataBase) -> None:
        self.db = db

    def update_db_table(self, table_name: str, table: pd.core.frame.DataFrame) -> None:
        self.db.tables[table_name] = table


class PermutationStatistic:

    rounding_decimal = 10
    def score_round(self, x): return round(x, self.rounding_decimal)

    def __init__(self):
        self.observed = None
        self.shuffles = {}

    def set_data(self, score, entry=None):
        self.observed = {'score': score, 'entry': entry}

    def add_shuffle_data(self, score, entry=None, seed=None, key=None):
        if not bool(key):
            key = len(self.shuffles)
        self.shuffles[key] = {'score': self.score_round(score), 'entry': entry, 'seed': seed}

    def get_score_distribution(self):
        return [shuffle['score'] for shuffle in self.shuffles.values()]

    def distribution_p_threshold(self, percentile):
        return np.percentile(self.get_score_distribution(), percentile, interpolation='linear')

    def significant(self, percentile):
        return self.observed['score'] > self.distribution_p_threshold(percentile)

    def zscore(self, key=None):
        shuffle_scores = self.get_score_distribution()
        val = self.observed['score'] if not bool(key) else self.shuffles[key['score']]
        return (val - np.mean(shuffle_scores) / np.std(shuffle_scores)) if np.std(shuffle_scores) > 0 else 0


class ObjectTimeseries:

    def __init__(self, timebin_interval: TimebinInterval=TimebinInterval(0, 0, 0, 0), series=None):
        self.timebin_interval = timebin_interval
        self.series = series

    def set_series(self, series):
        self.series = series

    def crop_timeseries(self, t_start, t_end):
        ''' t_start represents onset of first window and t_end represents offset of last window '''
        crop_start_ind = self.timebin_interval.split_to_bins_offset().index(t_start)
        crop_end_ind = self.timebin_interval.split_to_bins_offset().index(t_end)
        crop_slice = slice(crop_start_ind, crop_end_ind + 1)
        self.series = self.series[crop_slice]
        self.timebin_interval = self.timebin_interval.sub_interval(t_start - self.timebin_interval.timebin, t_end)

    def add_series_entry(self, entry, pos: int=-1):
        ''' -1 for add to end or 0 for start '''
        self.timebin_interval = copy(self.timebin_interval)
        if pos == -1:
            self.series.append(entry)
            self.timebin_interval.t_end += self.timebin_interval.timestep
        elif pos == 0:
            self.series.insert(0, entry)
            self.timebin_interval.t_start -= self.timebin_interval.timestep


class TableOperator:

    def __init__(self, df):
        self.df = df
        self.valid_assessed = None
        self.valid = None
        self.counts_min = None

    def assess_conditions(self, columns, levels, counts_thr, return_full=True):
        grouper = self.df.groupby(columns)
        group_len = grouper.apply(len)
        valid = set(levels) == set(grouper.groups.keys()) and (group_len >= counts_thr).all()
        counts_min = group_len.min()
        levels_groups = grouper.apply(lambda g: list(g.index)) if return_full else None
        self.valid_assessed = True
        self.valid = valid
        self.counts_min = counts_min
        return valid, counts_min, levels_groups

    def sample_conditions(self, columns, levels, counts_thr, max_sample=False, replace=False, random_state=None):
        if not self.valid_assessed:
            self.assess_conditions(columns, levels, counts_thr, return_full=False)
        if self.valid:
            n_samples = self.counts_min if max_sample else counts_thr
            return self.df.groupby(columns).sample(n_samples, random_state=random_state, replace=replace)
        else:
            return pd.Data.DataFrame([], columns=self.df.columns)


class BehavioralUnitFiringRate(Analysis):

    def __init__(self, db: DataBase, version: dict, unit_id: Union[int, tuple]=None) -> None:

        super(BehavioralUnitFiringRate, self).__init__(db)
        # reqs: units
        self.version = version
        self.data = None
        self.ind = None
        if type(unit_id) is int:
            self.load_from_iloc(unit_id)
        elif type(unit_id) is tuple:
            self.load_from_ind(unit_id)

    # TODO: add methods for constructing inspired from behunit_FR

    # IO
    def load_from_iloc(self, iloc: int) -> None:
        self.load_from_ind(self.db.tables['units'].iloc[iloc].name)

    def load_from_ind(self, ind: tuple) -> None:
        self.ind = ind
        self.data = MetaData().np_loader(self.get_data_filename(ind))

    def get_data_filename(self, ind: tuple) -> str:
        version_fr = self.version['fr']
        v_fr_params = anova_version_fr_params(version_fr)
        return MetaData().proc_dest_path(path.join('BehavioralUnits', 'FiringRates',
                                                   behunit_params_str(version_fr,
                                                                      v_fr_params['timebin'], v_fr_params['timestep'],
                                                                      v_fr_params['t_start'], v_fr_params['t_end'])),
                                         'behunit_FR_{0:s}_chan{1:03d}_unit{2:03d}.pkl'.format(*ind))


class PopulationBehavioralTimeseries:

    base_columns = ['Unit', 'Event', 'Condition', 'Instance', 'Timeseries']

    def __init__(self, condition_labels=None, timebin_interval=TimebinInterval(0, 0, 0, 0)):
        self.df = pd.DataFrame(columns=self.base_columns)
        self.condition_labels = [] if not bool(condition_labels) else condition_labels
        self.timebin_interval = timebin_interval

    def init_with_df(self, df):
        pbt = PopulationBehavioralTimeseries(condition_labels=self.condition_labels, timebin_interval=self.timebin_interval)
        pbt.df = df
        return pbt

    def get_unit_inds(self):
        return list(self.df['Unit'].unique())

    def get_condition_levels_1d(self):
        return list(self.df['Condition'].unique())

    def get_condition_levels_nd(self):
        cond_1d_list = self.get_condition_levels_1d()
        return [tuple(np.unique(t)) for t in zip(*cond_1d_list)]

    def get_condition_slice(self, condition_list):
        return self.df.loc[self.df['Condition'].isin(condition_list)]

    def get_unit_slice(self, unit_list):
        return self.df.loc[self.df['Unit'].isin(unit_list)]

    def add_data_row(self, unit_ind, event_ind, condition, instance, timeseries, index=None):
        row_data_df = pd.DataFrame.from_dict({index: [unit_ind, event_ind, condition, instance, timeseries]},
                                             columns=self.base_columns,
                                             orient='index')
        ignore_index = index is None
        self.df.append(row_data_df, ignore_index=ignore_index)

    def add_data_rows_from_list(self, data_list):
        row_data_df = pd.DataFrame(data_list, columns=self.base_columns)
        self.df = self.df.append(row_data_df, ignore_index=False)

    def func_instances(self, group_columns, func):
        grouper = self.df.groupby(group_columns)
        grouper_index = grouper.groups.keys()
        timeseries = grouper['Timeseries'].aggregate(func)
        events = grouper['Event'].aggregate(list)
        df_average = pd.DataFrame({'Timeseries': timeseries,
                                   'Event': events,
                                   'Instance': [0 for _ in range(len(grouper_index))]},
                                  index=grouper_index)
        df_average.index.set_names(group_columns, inplace=True)
        df_average.reset_index(drop=False, inplace=True)
        return self.init_with_df(df_average)

    def average_instances(self, group_columns):
        return self.func_instances(group_columns, lambda x: list(np.mean(np.array(list(x)), axis=0)))

    def std_instances(self, group_columns):
        return self.func_instances(group_columns, lambda x: list(np.std(np.array(list(x)), axis=0)))

    def se_instances(self, group_columns):
        return self.func_instances(group_columns, lambda x: list(np.std(np.array(list(x)), axis=0) / np.sqrt(len(x))))

    def crop_timeseries(self, t_start, t_end):
        crop_start_ind = self.timebin_interval.split_to_bins_offset().index(t_start)
        crop_end_ind = self.timebin_interval.split_to_bins_offset().index(t_end)
        crop_slice = slice(crop_start_ind, crop_end_ind + 1)
        self.df['Timeseries'] = self.df['Timeseries'].apply(lambda x: x[crop_slice])
        self.timebin_interval = self.timebin_interval.sub_interval(t_start - self.timebin_interval.timebin, t_end)

    def smooth_timeseries(self, smoother):
        return self.df['Timeseries'].apply(smoother.smoothen)

    def smooth_df(self, smoother):
        smooth_timeseries = self.smooth_timeseries(smoother)
        df_smooth = self.df.copy()
        df_smooth['Timeseries'] = smooth_timeseries
        return df_smooth

    def unfold_conditions_df(self):
        return pd.DataFrame(self.df['Condition'].tolist(), index=self.df.index, columns=self.condition_labels)

    def to_PCA_array(self):
        df_base = self.df
        grouper = df_base.groupby(['Unit', 'Condition', 'Instance'])
        num_units = df_base['Unit'].nunique()
        num_conditions = df_base['Condition'].nunique()
        num_instances = df_base.groupby(['Unit', 'Condition']).size().unique()[0]
        num_timebins = self.timebin_interval.num_of_bins()
        # order indices by groupping order unit, condition, instance
        index_list = [index for _, group in grouper for index in group.index]
        # create pre_pca unit, condition, instance, timebin array
        X_preshape = (num_units, num_conditions, num_instances, num_timebins)
        X_pre = np.array(df_base.iloc[index_list]['Timeseries'].tolist()).reshape(X_preshape)
        # convert to pca array format
        X_shape = (num_conditions * num_instances * num_timebins, num_units)
        X = X_pre.transpose(1, 2, 3, 0).reshape(X_shape)
        # unit, condition, instance information for reconstruction
        records = df_base.iloc[index_list][['Unit', 'Condition', 'Instance']].to_records()
        return X, records

    def to_PCA_scale_array(self):
        X, records = self.to_PCA_array()
        X_scale = StandardScaler().fit_transform(X.transpose()).transpose()
        return (X, X_scale), records

    def to_dPCA_mean_array(self):
        # init params
        pbt_mean = self.average_instances(['Unit', 'Condition'])
        df_base = pd.concat([pbt_mean.df, pbt_mean.unfold_conditions_df()], axis=1)
        grouper = df_base.groupby(['Unit'] + pbt_mean.condition_labels + ['Instance'])
        num_units = df_base['Unit'].nunique()
        num_conditions = list(df_base[pbt_mean.condition_labels].nunique())
        num_instances = df_base.groupby(['Unit', 'Condition']).size().unique()[0]
        num_timebins = pbt_mean.timebin_interval.num_of_bins()
        # order indices by groupping order unit, condition, instance
        index_list = [index for _, group in grouper for index in group.index]
        # create pre_pca unit, condition, instance, timebin array
        X_preshape = tuple([num_units] + num_conditions + [num_instances, num_timebins])
        X_pre = np.array(df_base.iloc[index_list]['Timeseries'].tolist()).reshape(X_preshape)
        # convert to pca array format
        X = X_pre.squeeze()
        # unit, condition, instance information for reconstruction
        records = df_base.iloc[index_list][['Unit', 'Condition', 'Instance']].to_records()
        return X, records

    def to_dPCA_mean_demean_array(self):
        X_mean, records = self.to_dPCA_mean_array()
        mean_shape = X_mean.shape
        X__mean_2d = X_mean.reshape((mean_shape[0], -1))
        X_mean_demean = StandardScaler(with_std=False).fit_transform(X__mean_2d.transpose()).transpose().reshape(mean_shape)
        return (X_mean, X_mean_demean), records

    def to_dPCA_trial_array(self):
        # init params
        df_base = pd.concat([self.df, self.unfold_conditions_df()], axis=1)
        grouper = df_base.groupby(['Unit'] + self.condition_labels + ['Instance'])
        num_units = df_base['Unit'].nunique()
        num_conditions = list(df_base[self.condition_labels].nunique())
        num_instances = df_base.groupby(['Unit', 'Condition']).size().unique()[0]
        num_timebins = self.timebin_interval.num_of_bins()
        # order indices by groupping order unit, condition, instance
        index_list = [index for _, group in grouper for index in group.index]
        # create pre_pca unit, condition, instance, timebin array
        X_preshape = tuple([num_units] + num_conditions + [num_instances, num_timebins])
        X_pre = np.array(df_base.iloc[index_list]['Timeseries'].tolist()).reshape(X_preshape)
        # convert to pca array format
        X_transpose_order = tuple([len(num_conditions) + 1, 0] + list(range(1, len(num_conditions) + 1)) + [len(num_conditions) + 2])
        X = X_pre.transpose(X_transpose_order)
        # unit, condition, instance information for reconstruction
        records = df_base.iloc[index_list][['Unit', 'Condition', 'Instance']].to_records()
        return X, records


class FactorBehavioralTimeseries(PopulationBehavioralTimeseries):

    base_columns = ['Factor', 'Condition', 'Instance', 'Timeseries']

    def __init__(self, df, condition_labels, timebin_interval):
        self.df = df
        self.condition_labels = condition_labels
        self.timebin_interval = timebin_interval

    def init_with_df(self, df):
        pbt = PopulationBehavioralTimeseries(condition_labels=self.condition_labels,
                                             timebin_interval=self.timebin_interval)
        pbt.df = df
        return pbt

    def unfold_conditions_df(self):
        return pd.DataFrame(self.df['Condition'].tolist(), index=self.df.index, columns=self.condition_labels)

    def average_instances(self, group_columns):
        grouper = self.df.groupby(group_columns)
        grouper_index = grouper.groups.keys()
        timeseries = grouper['Timeseries'].aggregate(lambda x: list(np.mean(np.array(list(x)), axis=0)))
        df_average = pd.DataFrame({'Timeseries': timeseries,
                                   'Instance': [0 for _ in range(len(grouper_index))]},
                                  index=grouper_index)
        df_average.index.set_names(group_columns, inplace=True)
        df_average.reset_index(drop=False, inplace=True)
        return FactorBehavioralTimeseries(df_average, self.condition_labels, self.timebin_interval)

    def init_with_df(self): pass
    def get_unit_inds(self): pass
    def get_unit_slice(self): pass
    def add_data_row(self): pass
    def add_data_rows_from_list(self): pass
    def to_PCA_array(self): pass
    def to_dPCA_mean_array(self): pass
    def to_dPCA_trial_array(self): pass


class DemixedPrincipalComponent(Analysis):

    def __init__(self, db: DataBase, version: dict) -> None:
        # reqs: trials, units, events, conditions
        self.db = db
        self.version = version
        self.pbt = None
        self.fbt = None
        self.dpca_obj = None

    def assess_unit_events(self, unit_id: Union[int, tuple]) -> list:

        bufr = BehavioralUnitFiringRate(self.db, {'fr': self.version['fr']}, unit_id)

        # functions for evaluating number of events for each condition for unit thresholded by mean firing rate
        def group_inds(group): return group.index
        def valid_timeseries_fr(timeseries): return np.mean(timeseries) < float(self.version['fr_thr'])
        def valid_events(inds): return [ind for ind in inds if valid_timeseries_fr(bufr.data[ind]) and self.db.tables['trials'].loc[ind[0:2]].StopCondition == 1]

        # params
        v_factor_params = factor_generate_conditions(self.version['factor'])
        condition_columns = v_factor_params['condition_columns']
        condition_list = v_factor_params['condition_list']

        # group unit's events by condition
        condition_grouper = self.db.tables['events_conditions'].loc[list(bufr.data.keys())].groupby(condition_columns)
        # if unit has events for all conditions and valid events are more than counts threshold return list of events
        events_inds = []
        for condition in condition_list:
            if condition in condition_grouper.groups.keys() and len(valid_events(group_inds(condition_grouper.get_group(condition)))) >= int(self.version['counts_thr']):
                events_inds.extend(valid_events(group_inds(condition_grouper.get_group(condition))))
            else:
                events_inds = []
                break
        return events_inds

    def exec(self) -> tuple:

        if not bool(self.pbt):
            self.pbt = self.db.md.np_loader(self.get_exec_filename('pbt'))

        X_trial, records_trial = self.pbt.to_dPCA_trial_array()
        (X_pre, X), records = self.pbt.to_dPCA_mean_demean_array()
        labels = factor_dpca_labels_mapping(self.version['factor'])
        join = factor_dpca_join_mapping(labels)
        self.dpca_obj = dPCA.dPCA(labels=labels, join=join, regularizer='auto')
        self.dpca_obj.n_trials = 5
        self.dpca_obj.protect = ['t']
        X_fit = self.dpca_obj.fit_transform(X, X_trial)
        self.dpca_obj.join = self.dpca_obj._join  # TODO: figure out why this throws bug if not added here

        return self.dpca_obj, X_fit, X, X_trial, records, records_trial

    # IO
    def get_wrangle_filename(self) -> str:
        version_fr = self.version['fr']
        version_factor = self.version['factor']
        counts_thr = int(self.version['counts_thr'])
        v_fr_params = anova_version_fr_params(version_fr)
        return MetaData().proc_dest_path(path.join('BehavioralUnits', 'DemixedPCA', version_factor,
                                                   behunit_params_str(version_fr,
                                                                      v_fr_params['timebin'], v_fr_params['timestep'],
                                                                      v_fr_params['t_start'], v_fr_params['t_end']),
                                                   '{0:03d}'.format(counts_thr), 'wrangle'),
                                         'valid_units.pkl')

    def get_filter_filename(self) -> str:
        version_fr = self.version['fr']
        version_factor = self.version['factor']
        counts_thr = int(self.version['counts_thr'])
        v_fr_params = anova_version_fr_params(version_fr)
        area_list_str = list_param_to_str(self.version['area_list'])
        subject_str = list_param_to_str(self.version['subject'])
        return MetaData().proc_dest_path(path.join('BehavioralUnits', 'DemixedPCA', version_factor,
                                                   behunit_params_str(version_fr,
                                                                      v_fr_params['timebin'], v_fr_params['timestep'],
                                                                      v_fr_params['t_start'], v_fr_params['t_end']),
                                                   '{0:03d}'.format(counts_thr), 'filter',
                                                   '_'.join([area_list_str, subject_str])),
                                         'filter.pkl')

    def get_exec_filename(self, filename) -> str:
        # 'pbt', 'fbt', 'dpca_obj', 'X_fit', 'X_tuple', 'significance'
        version_fr = self.version['fr']
        version_factor = self.version['factor']
        counts_thr = int(self.version['counts_thr'])
        v_fr_params = anova_version_fr_params(version_fr)
        area_list_str = list_param_to_str(self.version['area_list'])
        subject_str = list_param_to_str(self.version['subject'])
        return MetaData().proc_dest_path(path.join('BehavioralUnits', 'DemixedPCA', version_factor,
                                                   behunit_params_str(version_fr,
                                                                      v_fr_params['timebin'], v_fr_params['timestep'],
                                                                      v_fr_params['t_start'], v_fr_params['t_end']),
                                                   '{0:03d}'.format(counts_thr), 'exec',
                                                   '_'.join([area_list_str, subject_str]),
                                                   self.version['area'],
                                                   '_'.join([self.version['mode'], '{0:03d}'.format(int(self.version['mode_seed']))])),
                                         '{0:s}.pkl'.format(filename))


class ClassificationAnalysis(Analysis):

    def __init__(self, db: DataBase, version: dict) -> None:
        # reqs: trials, units, events, conditions
        self.db = db
        self.version = version

    def assess_unit_events(self, unit_id: Union[int, tuple]) -> list:

        bufr = BehavioralUnitFiringRate(self.db, {'fr': self.version['fr']}, unit_id)

        # params
        v_class_params = classification_version_class(self.version['class'])
        v_balance_params = classification_version_balance(self.version['balance'])
        condition_columns = v_class_params['condition_columns'] + v_balance_params['condition_columns']
        condition_list = list(product(*v_class_params['condition_list'], *v_balance_params['condition_list']))

        # filter only events of correct trials, observed by unit and belonging on condition_list
        events_conditions_db = self.db.tables['events_conditions']
        events_conditions_unit = events_conditions_db.loc[bufr.data.keys()]
        events_conditions = events_conditions_unit.loc[events_conditions_unit['StopCondition'].eq(1) & \
                                                       [el in condition_list for el in combine_columns(events_conditions_unit, condition_columns)]]
        ect_operator = TableOperator(events_conditions)
        # assess events, conditions, trials for counterbalanced number of condition list occuring in condition columns
        ect_assessment = ect_operator.assess_conditions(condition_columns, condition_list, int(self.version['counts_thr']), return_full=True)
        # if valid unit
        if ect_assessment[0]:
            # flatten df with list of events into a single list
            return list(chain(*list(ect_assessment[2])))
        else:
            return []

    # IO
    def get_wrangle_filename(self) -> str:
        version_fr = self.version['fr']
        version_factor = self.version['class']
        version_balance = self.version['balance']
        counts_thr = int(self.version['counts_thr'])
        v_fr_params = anova_version_fr_params(version_fr)
        return MetaData().proc_dest_path(path.join('BehavioralUnits', 'Classification', version_factor, version_balance,
                                                   behunit_params_str(version_fr,
                                                                      v_fr_params['timebin'], v_fr_params['timestep'],
                                                                      v_fr_params['t_start'], v_fr_params['t_end']),
                                                   '{0:03d}'.format(counts_thr), 'wrangle'),
                                         'valid_units.pkl')


# ### Misc ### #

def list_param_to_list(list_param):
    return sorted(list_param.split('_'))


def list_param_to_str(list_param):
    return ''.join(list_param_to_list(list_param))
