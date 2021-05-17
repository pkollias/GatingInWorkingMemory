from __future__ import annotations
from typing import Union
from rec_format import *
from versioning import *
from rec_db import *
from rec_stats import *
from dPCA import dPCA
from sklearn.preprocessing import StandardScaler


class Analysis:

    def __init__(self, db: DataBase) -> Analysis:
        self.db = db

    def update_db_table(self, table_name: str, table: pd.core.frame.DataFrame) -> None:
        self.db.tables[table_name] = table


class PermutationStatistic:

    rounding_decimal = 10
    score_round = lambda self, x: round(x, self.rounding_decimal)

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


class BehavioralUnitFR(Analysis):

    def __init__(self, db: DataBase, version: dict, unit_id: Union[int, tuple]=None) -> BehavioralUnitFR:
        # reqs: units
        self.db = db
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


class DemixedPrincipalComponent(Analysis):

    def __init__(self, db: DataBase, version: dict) -> DemixedPrincipalComponent:
        # reqs: trials, units, events, conditions
        self.db = db
        self.version = version
        self.pbt = None
        self.fbt = None

    def assess_unit_events(self, unit_id: Union[int, tuple]) -> list:

        bufr = BehavioralUnitFR(self.db, {'fr': self.version['fr']}, unit_id)

        # functions for evaluating number of events for each condition for unit thresholded by mean firing rate
        group_inds = lambda group: group.index
        valid_timeseries_fr = lambda timeseries: np.mean(timeseries) < float(self.version['fr_thr'])
        valid_events = lambda inds: [ind for ind in inds
                                     if valid_timeseries_fr(bufr.data[ind]) and
                                        self.db.tables['trials'].loc[ind[0:2]].StopCondition == 1]

        #params
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

    def pbt_from_behavioral_units(self, behavioral_units: pd.core.frame.DataFrame) -> PopulationBehavioralTimeseries:

        # params
        v_factor_params = factor_generate_conditions(self.version['factor'])
        condition_columns = v_factor_params['condition_columns']
        v_fr_params = anova_version_fr_params(self.version['fr'])
        t_start = v_fr_params['t_start']
        t_end = v_fr_params['t_end']
        timebin = v_fr_params['timebin']
        timestep = v_fr_params['timestep']

        # init
        timebin_interval = TimebinInterval(timebin, timestep, t_start, t_end)

        data_list = []
        # find units of behavioral units
        bu_grouper = behavioral_units.groupby(self.db.md.preproc_imports['units']['index'])
        # for every unit
        for unit_ind, unit_group in bu_grouper:

            # load firing rate data
            bufr = BehavioralUnitFR(self.db, {'fr': self.version['fr']}, unit_ind)

            # group by condition
            bu_cond_grouper = unit_group.groupby(condition_columns)
            # for every condition
            for condition, cond_group in bu_cond_grouper:
                # for every event
                for ii, (_, event_entry) in enumerate(cond_group.iterrows()):
                    event_ind = tuple(event_entry[self.db.md.preproc_imports['events']['index']])
                    condition = tuple(event_entry[condition_columns])
                    timeseries = bufr.data[event_ind]
                    data_entry = [unit_ind, event_ind, condition, ii, timeseries]
                    # append to data_list
                    data_list.append(data_entry)

        # create pbt
        pbt = PopulationBehavioralTimeseries(condition_columns, timebin_interval)
        pbt.add_data_rows_from_list(data_list)

        return pbt

    def exec(self, pbt: PopulationBehavioralTimeseries) -> tuple:

        X_trial, records_trial = pbt.to_dPCA_trial_array()
        X_pre, records = pbt.to_dPCA_mean_array()
        X = self.demean_mean_X(X_pre)
        labels = factor_dpca_labels_mapping(self.version['factor'])
        join = factor_dpca_join_mapping(labels)
        dpca_obj = dPCA.dPCA(labels=labels, join=join, regularizer='auto')
        dpca_obj.n_trials = 5
        dpca_obj.protect = ['t']
        X_fit = dpca_obj.fit_transform(X, X_trial)
        dpca_obj.join = dpca_obj._join # TODO: figure out why this throws bug if not added here

        return dpca_obj, X_fit, X, X_trial, records, records_trial

    def demean_mean_X(self, X_pre: np.ndarray) -> np.ndarray:

        X_shape = X_pre.shape
        X_pre_2d = X_pre.reshape((X_shape[0], -1))
        X = StandardScaler(with_std=False).fit_transform(X_pre_2d.transpose()).transpose().reshape(X_shape)
        return X

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


# ### Misc ### #

def list_param_to_list(list_param):
    return sorted(list_param.split('_'))

def list_param_to_str(list_param):
    return ''.join(list_param_to_list(list_param))
