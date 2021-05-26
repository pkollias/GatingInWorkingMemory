from __future__ import annotations
from metadata import *


class DataBase:

    def __init__(self, base: list = None) -> DataBase:
        self.md = MetaData()
        self.tables = self.md.db_base_loader(base)
        self.tables_src = self.tables
        if all([table in self.tables.keys() for table in ['events', 'conditions']]):
            self.tables['events_conditions'] = self.generate_events_conditions()

    def reset_tables(self) -> None:
        self.tables = self.tables_src

    def get_table_index(self, table_name: str) -> list:
        if table_name in self.md.preproc_imports.keys():
            return self.md.preproc_imports[table_name]['index']
        elif table_name in self.md.proc_imports.keys():
            return self.md.proc_imports[table_name]['index']

    def filter_db_by_table_inds(self, table_name: str, inds: list) -> None:
        self.tables[table_name] = self.tables_src[table_name].loc[inds]
        self.tighten_db()

    # TODO: implement later, backwards and forwards tightening based on entries left on tables
    def tighten_db(self):
        pass

    # units

    def unit_ind_from_iloc(self, ii: int) -> tuple:
        units = self.tables['units']
        return units.iloc[ii].name

    def unit_is_multiunit(self, ind: tuple) -> bool:
        units = self.tables['units']
        return units.loc[ind].UnitNum == 0 or units.loc[ind].RatingCode == 7

    # events, conditions

    def generate_events_conditions(self, src: bool = False) -> pd.core.frame.DataFrame:
        events_index = self.md.preproc_imports['events']['index']
        events = self.tables_src['events'].reset_index(drop=True) if src else self.tables['events'].reset_index(drop=True)
        conditions = self.tables_src['conditions'] if src else self.tables['conditions']
        events_conditions = pd.merge(events, conditions, on=events_index)
        events_conditions.set_index(events_index, drop=False, inplace=True)
        return events_conditions

    # events, conditions, trials

    def merge_events_conditions_trials(self, src: bool = False):
        trials_index = self.md.preproc_imports['trials']['index']
        events_index = self.md.preproc_imports['events']['index']
        self.tables['events_conditions'] = pd.merge(self.tables['events_conditions'].reset_index(drop=True),
                                                    self.tables['trials'].reset_index(drop=True),
                                                    on=trials_index).set_index(events_index)
    # activity

    def timestamp_interval_within_activity(self, start, end):
        activity = self.tables['activity']
        return (activity['SegmentStart'].le(start) & activity['SegmentEnd'].gt(start)).any() and \
               (activity['SegmentStart'].le(end) & activity['SegmentEnd'].gt(end)).any()



# ### Misc ### #


def combine_columns(table, columns):
    return list(zip(*[table[col] for col in columns]))


def timestamp_interval_within_activity(interval_start, interval_end, activity_collection):

    return (activity_collection['SegmentStart'].le(interval_start) &
            activity_collection['SegmentEnd'].gt(interval_start)).any() \
           and \
           (activity_collection['SegmentStart'].le(interval_end) &
            activity_collection['SegmentEnd'].gt(interval_end)).any()
