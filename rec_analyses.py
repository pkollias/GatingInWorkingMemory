from __future__ import annotations
from typing import Union
from rec import TimebinInterval
from rec_db import *
from rec_format import *
from versioning import *
from operator import itemgetter
from dPCA import dPCA
from copy import copy
from itertools import chain
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split


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

    def get_score_distribution(self) -> list:
        return [shuffle['score'] for shuffle in self.shuffles.values()]

    def distribution_p_threshold(self, percentile: float) -> float:
        return np.percentile(self.get_score_distribution(), percentile, interpolation='linear')

    def significant(self, percentile: float) -> bool:
        return self.observed['score'] > self.distribution_p_threshold(percentile)

    def zscore(self, key=None) -> float:
        shuffle_scores = self.get_score_distribution()
        val = self.observed['score'] if not bool(key) else self.shuffles[key['score']]
        return (val - np.mean(shuffle_scores) / np.std(shuffle_scores)) if np.std(shuffle_scores) > 0 else 0


class ObjectTimeseries:

    def __init__(self, timebin_interval: TimebinInterval=TimebinInterval(0, 0, 0, 0), series=None):
        self.timebin_interval = timebin_interval
        self.series = series

    def set_series(self, series):
        self.series = series

    def crop_timeseries(self, t_start: int, t_end: int):
        """ t_start represents onset of first window and t_end represents offset of last window """
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

    def get_object_at_time(self, t):
        return self.series[self.timebin_interval.t_offset_to_ind(t)]


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


class Anova(Analysis):

    def __init__(self, db: DataBase, version: dict) -> None:
        self.db = db
        self.version = version
        self.physiology_dict = None

    def load_physiology_dict(self, cropped=True) -> None:
        cropped_fname = self.get_path_base('physiology_dict_cropped', self.get_summary_consolidate_stem(cluster_corrected=True))
        uncropped_fname = self.get_path_base('physiology_dict', self.get_summary_consolidate_stem(cluster_corrected=True))
        if path.exists(cropped_fname) and cropped:
            print(cropped_fname)
            self.physiology_dict = self.db.md.np_loader(cropped_fname)
        else:
            print(uncropped_fname)
            self.physiology_dict = self.db.md.np_loader(uncropped_fname)

    def get_all_units(self, selection: str, x_factor: str) -> list:
        if not bool(self.physiology_dict):
            self.load_physiology_dict()
        return [unit_ind for unit_ind, unit_entry in self.physiology_dict[selection].items()
                if unit_entry[x_factor]['valid']]

    def get_selective_units(self, selection: str, x_factor: str) -> list:
        return [unit_ind for unit_ind, unit_entry in self.physiology_dict[selection].items()
                if unit_entry[x_factor]['valid'] and bool(unit_entry[x_factor]['clusters'])]

    def get_unit_time_selectivity(self, selection: str, x_factor: str, timebin_interval: TimebinInterval, filtered_inds: list=None) -> np.ndarray:
        # init
        if not bool(self.physiology_dict):
            self.load_physiology_dict()
        if not filtered_inds:
            filtered_inds = self.physiology_dict[selection].keys()

        tb_keys = timebin_interval.split_to_bins_onset()
        selectivity_array = np.zeros((len(filtered_inds), len(tb_keys)))
        for unit_ii, unit_ind, in enumerate(filtered_inds):
            physiology_entry = self.physiology_dict[selection][unit_ind]
            if not physiology_entry[x_factor]['valid']:
                selectivity_array[unit_ii, :] = np.nan
            else:
                selectivity_array[unit_ii, [tb_keys.index(el) for cl in physiology_entry[x_factor]['clusters'] for el in cl]] = 1

        return selectivity_array

    # IO
    def get_path_base(self, filename, stem='') -> str:
        version_fr = self.version['fr']
        version_aov = self.version['aov']
        v_fr_params = version_fr_params(version_fr)
        if type(stem) == str:
            stem_path = stem
        elif type(stem) == tuple:
            stem_path = path.join(*stem)
        return MetaData().proc_dest_path(path.join('BehavioralUnits', 'Anova', version_aov,
                                                   behunit_params_str(self.version['fr'], v_fr_params['timebin'], v_fr_params['timestep'],
                                                                      v_fr_params['t_start'], v_fr_params['t_end']),
                                                   stem_path),
                                         '{0:s}.pkl'.format(filename))

    def get_summary_consolidate_stem(self, cluster_corrected: bool) -> Union[str, tuple]:
        return 'consolidate' if cluster_corrected else 'summarize'


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
        v_fr_params = version_fr_params(version_fr)
        return MetaData().proc_dest_path(path.join('BehavioralUnits', 'FiringRates',
                                                   behunit_params_str(version_fr,
                                                                      v_fr_params['timebin'], v_fr_params['timestep'],
                                                                      v_fr_params['t_start'], v_fr_params['t_end'])),
                                         'behunit_FR_{0:s}_chan{1:03d}_unit{2:03d}.pkl'.format(*ind))


class PopulationBehavioralTimeseries:

    base_columns = ['Unit', 'Event', 'Condition', 'Instance', 'Timeseries']

    def __init__(self, condition_labels=None, timebin_interval=TimebinInterval(0, 0, 0, 0), df=None):
        if df is not None:
            self.df = df
        else:
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

    def to_behavioral_units(self):
        md = MetaData()
        unzipped_df = pd.concat([unzip_columns(self.df, 'Unit', md.preproc_imports['units']['index']),
                                 unzip_columns(self.df, 'Event', md.preproc_imports['events']['index']),
                                 unzip_columns(self.df, 'Condition', self.condition_labels)], axis=1)
        return unzipped_df.loc[:, ~unzipped_df.columns.duplicated()]

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

    def to_pseudosession_firing_rate(self, pseudosession_inds_array, t_ind):

        df_ordered_timeseries = self.df.loc[pseudosession_inds_array.flat]['Timeseries']
        df_t_fr = df_ordered_timeseries.apply(lambda timeseries: timeseries[t_ind])
        X = df_t_fr.to_numpy().reshape(pseudosession_inds_array.shape)
        return X


class PseudoPopulationBehavioralTimeseries(PopulationBehavioralTimeseries):

    base_columns = ['Unit', 'Unit_Code', 'Event', 'Condition', 'Instance', 'Timeseries']

    def __init__(self, condition_labels=None, timebin_interval=TimebinInterval(0, 0, 0, 0), df=None):
        if df is not None:
            self.df = df
        else:
            self.df = pd.DataFrame(columns=self.base_columns)
        self.condition_labels = [] if not bool(condition_labels) else condition_labels
        self.timebin_interval = timebin_interval

    def add_data_row(self, unit_ind, unit_code, event_ind, condition, instance, timeseries, index=None):
        row_data_df = pd.DataFrame.from_dict({index: [unit_ind, unit_code, event_ind, condition, instance, timeseries]},
                                             columns=self.base_columns,
                                             orient='index')
        ignore_index = index is None
        self.df.append(row_data_df, ignore_index=ignore_index)

    def to_behavioral_units(self):
        md = MetaData()
        unzipped_df = pd.concat([unzip_columns(self.df, 'Unit', md.preproc_imports['units']['index']),
                                 self.df['Unit_Code'],
                                 unzip_columns(self.df, 'Event', md.preproc_imports['events']['index']),
                                 self.df['Instance'],
                                 unzip_columns(self.df, 'Condition', self.condition_labels)], axis=1)
        return unzipped_df.loc[:, ~unzipped_df.columns.duplicated()]

    def to_PCA_array(self): pass
    def to_PCA_scale_array(self): pass
    def to_dPCA_mean_array(self): pass
    def to_dPCA_mean_demean_array(self): pass
    def to_dPCA_trial_array(self): pass


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
            if condition in condition_grouper.groups.keys() and len(condition_grouper.groups[condition]) \
                    and len(valid_events(group_inds(condition_grouper.get_group(condition)))) >= int(self.version['counts_thr']):
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
        self.dpca_obj.n_components = X_trial.shape[1]
        self.dpca_obj.protect = ['t']
        X_fit = self.dpca_obj.fit_transform(X, X_trial)
        self.dpca_obj.join = self.dpca_obj._join  # TODO: figure out why this throws bug if not added here

        return self.dpca_obj, X_fit, X, X_trial, records, records_trial

    # IO
    def get_path_base(self, filename, stem='') -> str:
        version_fr = self.version['fr']
        version_factor = self.version['factor']
        counts_thr = int(self.version['counts_thr'])
        v_fr_params = version_fr_params(version_fr)
        if type(stem) == str:
            stem_path = stem
        elif type(stem) == tuple:
            stem_path = path.join(*stem)
        return MetaData().proc_dest_path(path.join('BehavioralUnits', 'DemixedPCA', version_factor,
                                                   behunit_params_str(version_fr,
                                                                      v_fr_params['timebin'], v_fr_params['timestep'],
                                                                      v_fr_params['t_start'], v_fr_params['t_end']),
                                                   '{0:03d}'.format(counts_thr), stem_path),
                                         '{0:s}.pkl'.format(filename))

    def get_wrangle_stem(self) -> Union[str, tuple]:
        return 'wrangle'

    def get_filter_stem(self) -> Union[str, tuple]:
        area_list_str = list_param_to_str(self.version['area_list'])
        subject_str = list_param_to_str(self.version['subject'])
        return 'filter', '_'.join([area_list_str, subject_str])

    def get_exec_stem(self) -> Union[str, tuple]:
        # 'pbt', 'fbt', 'dpca_obj', 'X_fit', 'X_tuple', 'significance'
        area_list_str = list_param_to_str(self.version['area_list'])
        subject_str = list_param_to_str(self.version['subject'])
        return ('exec', '_'.join([area_list_str, subject_str]), self.version['area'],
                '_'.join([self.version['mode'], '{0:03d}'.format(int(self.version['mode_seed']))]))


class ClassificationAnalysis(Analysis):

    def __init__(self, db: DataBase, version: dict) -> None:
        # reqs: trials, units, events, conditions
        self.db = db
        self.version = version

    def get_assembly_condition_columns(self) -> list:

        v_class_params = classification_version_class(self.version['class'])
        v_balance_params = classification_version_balance(self.version['balance'])
        return v_class_params['condition_columns'] + v_balance_params['condition_columns']

    def get_assembly_condition_list(self) -> list:

        v_class_params = classification_version_class(self.version['class'])
        v_balance_params = classification_version_balance(self.version['balance'])
        return list(product(*v_class_params['condition_list'], *v_balance_params['condition_list']))

    def get_class_condition_columns(self) -> list:

        return classification_version_class(self.version['class'])['condition_columns']

    def get_class_condition_list(self) -> list:

        return list(product(*classification_version_class(self.version['class'])['condition_list']))

    def assess_unit_events(self, unit_id: Union[int, tuple]) -> list:

        bufr = BehavioralUnitFiringRate(self.db, {'fr': self.version['fr']}, unit_id)

        # params
        condition_columns = self.get_assembly_condition_columns()
        condition_list = self.get_assembly_condition_list()

        # filter only events of correct trials, observed by unit and belonging on condition_list
        events_conditions_db = self.db.tables['events_conditions']
        events_conditions_unit = events_conditions_db.loc[bufr.data.keys()]
        # get only correct trials and trials belonging in condition_list
        events_conditions = events_conditions_unit.loc[events_conditions_unit['StopCondition'].eq(1) & \
                                                       [el in condition_list for el in zip_columns(events_conditions_unit, condition_columns)]]
        ect_operator = TableOperator(events_conditions)
        # assess events, conditions, trials for counterbalanced number of condition list occuring in condition columns
        ect_assessment = ect_operator.assess_conditions(condition_columns, condition_list, int(self.version['counts_thr']), return_full=True)
        # if valid unit
        if ect_assessment[0]:
            # flatten df with list of events into a single list
            return list(chain(*list(ect_assessment[2])))
        else:
            return []

    def generate_pseudoarray_inds(self, pbt: PopulationBehavioralTimeseries, seed) -> tuple(np.ndarray, np.ndarray):

        # TODO: revisit - bad coding
        # important note: the shuffling happens at the level of the subgroups so internal structure and order is
        # preserved and not shuffled (e.g. order of output will be still S11, Gating -> S12, Gating, ...) so
        # class order can assumed to not be shuffled (therefore class order is still the same)

        md = MetaData()
        behavioral_units_index = md.proc_imports['behavioral_units']['index']
        units_index = md.preproc_imports['units']['index']

        # convert to df with beh_units columns and all other pbt.df columns unzipped
        # (session, channum, unitnum, trialnum, stageindex, + *condition + [unit_code, instance])
        behavioral_units = pbt.to_behavioral_units()
        pseudotrial_result = generate_pseudotrial_from_behavioral_units(behavioral_units,
                                                                        self.get_assembly_condition_columns(),
                                                                        self.get_class_condition_columns(),
                                                                        int(get_seed(hash((seed, '_'.join(self.version.values()))))))
        pseudotrial, unit_inds, unit_codes, class_codes, classes = pseudotrial_result

        # get C ordered list of behavioral units for X array
        def tuple_ind_to_behavioral_ind(tuple_ind): return (tuple_ind[0][0], tuple_ind[0][1], tuple_ind[0][2],  # unit_ind
                                                            tuple_ind[1],  # unit_code
                                                            tuple_ind[2][1], tuple_ind[2][2], tuple_ind[2][3])  # pseudo_event
        C_order_bu_inds = [tuple_ind_to_behavioral_ind(ind_tuple) for instance in pseudotrial for ind_tuple in zip(unit_inds, unit_codes, instance)]

        # TODO: create new index for behavioral pseudounits, amended with Unit_Code and Instance, correct later
        behavioral_pseudo_units_index = behavioral_units_index.copy()
        behavioral_pseudo_units_index.insert(len(units_index), 'Unit_Code')
        behavioral_pseudo_units_index.append('Instance')
        # transfer from C ordered list of behavioral unit indices to pbt int indices
        behavioral_unit_inds = zip_columns(pbt.to_behavioral_units(), behavioral_pseudo_units_index, 'Behavioral_PseudoUnit')
        C_order_pbt_iloc_inds = behavioral_unit_inds.reset_index(drop=False).set_index('Behavioral_PseudoUnit').loc[C_order_bu_inds]['index']
        pseudoarray_flat_inds = pbt.df.loc[C_order_pbt_iloc_inds].index.to_numpy()
        pseudoarray_inds = pseudoarray_flat_inds.reshape(len(classes), len(unit_inds))
        class_array = np.array(class_codes)

        return pseudoarray_inds, class_array

    def get_train_test_splits_list(self, pbt: PopulationBehavioralTimeseries):

        split = self.version['split']
        num_conditions_to_stratify = len(self.get_assembly_condition_list())
        instances_per_condition = int(self.version['counts_thr'])
        n_trials = num_conditions_to_stratify * instances_per_condition

        if split == 'StratifiedStim':

            # under design assumptions, X_inds_array holds latent (assembly) condition information
            # in order so we can tile stratify parameter
            condition_stratify = np.tile(np.arange(num_conditions_to_stratify),
                                         (instances_per_condition, 1)).transpose().reshape(-1)
            train_test = [train_test_split(np.arange(n_trials), train_size=2 / 3, stratify=condition_stratify, random_state=ii) for ii in
                          range(20)]

        elif split in ['OneStimTest', 'OneStimTrain', 'WithinGroupTransfer', 'AcrossGroupTransfer']:

            condition = unzip_columns(pbt.df, 'Condition', pbt.condition_labels)
            generalize_column = 'StageStimSpecialized'
            split_to_test_stim_levels = dict([('OneStimTest', [['S11'], ['S12'], ['S21'], ['S22']]),
                                              ('OneStimTrain', [['S11', 'S12', 'S21'], ['S11', 'S12', 'S22'], ['S11', 'S21', 'S22'], ['S12', 'S21', 'S22']]),
                                              ('WithinGroupTransfer', [['S11', 'S21'], ['S11', 'S22'], ['S12', 'S21'], ['S12', 'S22']]),
                                              ('AcrossGroupTransfer', [['S11', 'S12'], ['S21', 'S22']])])
            train_test = []
            for test_stim_list in split_to_test_stim_levels[split]:
                pbt_test_index_bool = condition.iloc[:n_trials][generalize_column].isin(test_stim_list)
                train_inds = np.random.permutation(np.array(pbt_test_index_bool.loc[~pbt_test_index_bool].index))
                test_inds = np.random.permutation(np.array(pbt_test_index_bool.loc[pbt_test_index_bool].index))
                train_test.append([train_inds, test_inds])

        elif split in ['StratifiedBalanceSplit', 'StratifiedBalanceSplit_StimHalf']:

            # under design assumptions, X_inds_array holds latent (assembly) condition information
            # in order so we can tile stratify parameter
            class_columns = classification_version_class(self.version['class'])['condition_columns']
            num_classes = len(np.unique(pbt.unfold_conditions_df()[class_columns].to_numpy()))
            assembly_balance_columns = classification_version_balance(self.version['balance'])['condition_columns']
            assembly_balance_list = classification_version_balance(self.version['balance'])['condition_list']
            split_ignore_columns = classification_version_balance(self.version['balance'])['split_ignore_columns']
            split_ignore_columns_index = list(map(assembly_balance_columns.index, split_ignore_columns))
            split_ignore_columns_list = [assembly_balance_list[ind] for ind in split_ignore_columns_index]
            num_values_to_ignore = list(map(len, split_ignore_columns_list))
            # assumption that split_columns have to be last in sequence within assembly balance columns
            balance_columns = [col for col in assembly_balance_columns if col not in split_ignore_columns]

            sub_pbt_df = zip_columns(pbt.unfold_conditions_df(), balance_columns, 'Balance_Condition') \
                if len(balance_columns) > 1 \
                else pbt.unfold_conditions_df()[balance_columns]
            conditions_to_split = np.unique(sub_pbt_df.to_numpy())
            if split in ['StratifiedBalanceSplit']:
                num_conditions_to_split = len(conditions_to_split)
                group_splits = np.arange(n_trials).reshape(num_classes, num_conditions_to_split, -1).transpose((1, 0, 2)).reshape(num_conditions_to_split, -1)
                condition_stratify = np.tile(np.arange(num_classes * np.product(num_values_to_ignore)),
                                             (instances_per_condition, 1)).transpose().reshape(-1)
            elif split in ['StratifiedBalanceSplit_StimHalf']:
                conditions_to_split = list(product(conditions_to_split, [1, 2]))
                num_conditions_to_split = len(conditions_to_split)
                group_splits = np.arange(n_trials).reshape(num_classes, int(num_conditions_to_split / 2), -1).transpose((1, 0, 2)).reshape(num_conditions_to_split, -1)
                condition_stratify = np.tile(np.arange(int((num_classes * np.product(num_values_to_ignore)) / 2)),
                                             (instances_per_condition, 1)).transpose().reshape(-1)

            train_test = {condition: [train_test_split(group_split, train_size=2 / 3,
                                                       stratify=condition_stratify, random_state=ii)
                                      for ii in range(20)]
                          for cond_ind, (condition, group_split)
                          in enumerate(zip(conditions_to_split, group_splits))}

        return train_test

    # IO
    def get_path_base(self, filename, stem='', cross=None) -> str:
        version_fr = self.version['fr']
        version_class = self.version['class'] if not cross else ''.join(sorted(self.version['class_list'].split('_'))) # TODO Correct !!! to join the split version
        version_balance = self.version['balance'] if not cross else ''.join(sorted(self.version['balance_list'].split('_'))) # TODO Correct !!! to join the split version
        counts_thr = int(self.version['counts_thr'])
        v_fr_params = version_fr_params(version_fr)
        if type(stem) == str:
            stem_path = stem
        elif type(stem) == tuple:
            stem_path = path.join(*stem)
        return MetaData().proc_dest_path(path.join('BehavioralUnits', 'Classification', version_class, version_balance,
                                                   behunit_params_str(version_fr,
                                                                      v_fr_params['timebin'], v_fr_params['timestep'],
                                                                      v_fr_params['t_start'], v_fr_params['t_end']),
                                                   '{0:03d}'.format(counts_thr), stem_path),
                                         '{0:s}.pkl'.format(filename))


    def get_wrangle_stem(self) -> Union[str, tuple]:
        return 'wrangle'

    def get_filter_stem(self) -> Union[str, tuple]:
        area_list_str = list_param_to_str(self.version['area_list'])
        subject_str = list_param_to_str(self.version['subject'])
        return 'filter', '_'.join([area_list_str, subject_str])
        return 'wrangle'

    def get_assemble_stem(self) -> Union[str, tuple]:
        area_list_str = list_param_to_str(self.version['area_list'])
        subject_str = list_param_to_str(self.version['subject'])
        return ('exec', '_'.join([area_list_str, subject_str]), self.version['area'],
                '_'.join([self.version['mode'], '{0:03d}'.format(int(self.version['mode_seed']))]))

    def get_train_stem(self) -> Union[str, tuple]:
        return (*self.get_assemble_stem(), 'train')

    def get_pseudo_session_stem(self) -> Union[str, tuple]:
        return (*self.get_train_stem(), 'pseudo_session_{0:03d}'.format(int(self.version['pseudo_seed'])))

    def get_train_test_stem(self) -> Union[str, tuple]:
        return (*self.get_pseudo_session_stem(), 'train_test', 'split_{0:s}'.format(self.version['split']))

    def get_filter_session_stem(self) -> Union[str, tuple]:
        area_list_str = list_param_to_str(self.version['area_list'])
        subject_str = list_param_to_str(self.version['subject'])
        return 'filter', '_'.join([area_list_str, subject_str, self.version['sess_ratio'], self.version['units_ratio']])

    def get_session_stem(self) -> Union[str, tuple]:
        return (*self.get_filter_session_stem(), 'session')

    def get_train_test_session_stem(self) -> Union[str, tuple]:
        return (*self.get_session_stem(), 'train_test')


# ### Misc ### #

def list_param_to_list(list_param: str) -> list:
    return sorted(list_param.split('_'))


def list_param_to_str(list_param: str) -> str:
    return ''.join(list_param_to_list(list_param))


def timebin_interval_from_version_fr(version_fr: str) -> TimebinInterval:
    v_fr_params = version_fr_params(version_fr)
    t_start = v_fr_params['t_start']
    t_end = v_fr_params['t_end']
    timebin = v_fr_params['timebin']
    timestep = v_fr_params['timestep']
    return TimebinInterval(timebin, timestep, t_start, t_end)
