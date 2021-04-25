from __future__ import annotations
from typing import Union
from rec_format import *
from rec import *
from rec_db import *


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

    class Params:
        def __init__(self, version_fr: str) -> BehavioralUnitFR.Params:
            self.version_fr = version_fr

    def __init__(self, db: DataBase, params: tuple, unit_id: Union[int, tuple]=None) -> BehavioralUnitFR:
        # reqs: units
        self.db = db
        self.params = BehavioralUnitFR.Params(*params) if type(params) is tuple else BehavioralUnitFR.Params(params)
        self.data = None
        self.ind = None
        if type(unit_id) is int:
            self.load_from_iloc(unit_id)
        elif type(unit_id) is tuple:
            self.load_from_ind(unit_id)

    # IO
    def load_from_iloc(self, iloc: int) -> None:
        self.load_from_ind(self.db.tables['units'].iloc[iloc].name)

    def load_from_ind(self, ind: int) -> None:
        self.ind = ind
        self.data = MetaData().np_loader(self.get_data_filename(ind))

    def get_data_filename(self, ind: tuple) -> str:
        version_fr = self.params.version_fr
        v_fr_params = anova_version_fr_params(version_fr)
        return MetaData().proc_dest_path(path.join('BehavioralUnits', 'FiringRates',
                                                   behunit_params_str(version_fr,
                                                                      v_fr_params['timebin'], v_fr_params['timestep'],
                                                                      v_fr_params['t_start'], v_fr_params['t_end'])),
                                         'behunit_FR_{0:s}_chan{1:03d}_unit{2:03d}.pkl'.format(*ind))


class DemixedPrincipalComponent(Analysis):

    class Params:
        def __init__(self, version_factor: str, version_fr: str, counts_thr: int = 20, fr_thr: float = 100) -> DemixedPrincipalComponent.Params:
            self.version_factor = version_factor
            self.version_fr = version_fr
            self.counts_thr = counts_thr
            self.fr_threshold = fr_thr

    def __init__(self, db: DataBase, params: tuple) -> DemixedPrincipalComponent:
        # reqs: units, events, conditions
        self.db = db
        self.params = DemixedPrincipalComponent.Params(*params) if type(params) is tuple else DemixedPrincipalComponent.Params(*params)
        self.pbt = None
        self.fbt = None

    # def construct_behavioral_units(self):
    #     beh_unit_inds_from_events_inds = lambda inds, unit_ind: [(*unit_ind, ind[1], ind[2]) for ind in inds]
    #     for unit_ind in units.index:
    #     behavioral_units_list = pd.DataFrame

    def assess_unit_events(self, unit_id: Union[int, tuple]) -> list:

        bufr = BehavioralUnitFR(self.db, (self.params.version_fr), unit_id)

        # functions for evaluating number of events for each condition for unit thresholded by mean firing rate
        group_inds = lambda group: group.index
        valid_timeseries_fr = lambda timeseries: np.mean(timeseries) < self.params.fr_threshold
        valid_events = lambda inds: [ind for ind in inds if valid_timeseries_fr(bufr.data[ind])]

        #params
        v_factor_params = factor_generate_conditions(self.params.version_factor)
        condition_columns = v_factor_params['condition_columns']
        condition_list = v_factor_params['condition_list']

        # group unit's events by condition
        condition_grouper = self.db.tables['events_conditions'].loc[list(bufr.data.keys())].groupby(condition_columns)
        # if unit has events for all conditions and valid events are more than counts threshold return list of events
        events_inds = []
        for condition in condition_list:
            if condition in condition_grouper.groups.keys() and len(valid_events(group_inds(condition_grouper.get_group(condition)))) >= self.params.counts_thr:
                events_inds.extend(valid_events(group_inds(condition_grouper.get_group(condition))))
            else:
                events_inds = []
                break
        return events_inds

    # IO
    def get_wrangle_filename(self) -> str:
        version_fr = self.params.version_fr
        version_factor = self.params.version_factor
        v_fr_params = anova_version_fr_params(version_fr)
        return MetaData().proc_dest_path(path.join('BehavioralUnits', 'DemixedPCA', version_factor,
                                                   behunit_params_str(version_fr,
                                                                      v_fr_params['timebin'], v_fr_params['timestep'],
                                                                      v_fr_params['t_start'], v_fr_params['t_end']),
                                                   'wrangle'),
                                         'valid_units.pkl')
