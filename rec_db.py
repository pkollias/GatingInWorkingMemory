from __future__ import annotations
from metadata import *
from rec_format import get_seed


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
                                                    on=trials_index, suffixes=('', '_trials')).set_index(events_index)
    # activity

    def timestamp_interval_within_activity(self, start, end):
        activity = self.tables['activity']
        return (activity['SegmentStart'].le(start) & activity['SegmentEnd'].gt(start)).any() and \
               (activity['SegmentStart'].le(end) & activity['SegmentEnd'].gt(end)).any()


# ### Misc ### #


def zip_columns(table, old_column_names, new_column=None):
    return pd.Series(list(zip(*[table[col] for col in old_column_names])), index=table.index, name=new_column)


def unzip_columns(table, old_column, new_column_names):
    return pd.DataFrame(list(map(list, table[old_column])), index=table.index, columns=new_column_names)


def timestamp_interval_within_activity(interval_start, interval_end, activity_collection):

    if np.isnan(interval_start) or np.isnan(interval_end):
        return False
    else:
        return (activity_collection['SegmentStart'].le(interval_start) &
                activity_collection['SegmentEnd'].gt(interval_start)).any() \
               and \
               (activity_collection['SegmentStart'].le(interval_end) &
                activity_collection['SegmentEnd'].gt(interval_end)).any()


def generate_pseudotrial_from_behavioral_units(behavioral_units, assembly_condition_columns, class_condition_columns, seed):

    '''groups behavioral units by pseudounit (unit_index + unit_code)
    within each pseudounit groupping groups by assembly_condition and shuffles order
    returns shuffled pseudounit, assembly_condition shuffles
    and information in preserved order of unit_inds, unit_codes, class_codes, classes'''

    md = MetaData()
    units_index = md.preproc_imports['units']['index']
    events_index = md.preproc_imports['events']['index']

    # util functions for shuffled pseudotrial generation
    def zip_events(df): return zip_columns(df, events_index + ['Instance'], 'PseudoEvent')
    # ensure that each group has its own seed
    def shuffle_group_rows(group): return group.groupby(assembly_condition_columns).sample(frac=1, random_state=get_seed((seed, group.name)))
    def apply_func(group): return list(zip_events(shuffle_group_rows(group)))

    # get transpose of pseudotrial event indices by shuffling events within condition and then stacking together
    unit_grouper = behavioral_units.groupby(units_index + ['Unit_Code'])
    pseudo_trials_t = list(unit_grouper.apply(apply_func))
    # transpose generated pseudotrial list
    pseudo_trials = list(map(list, zip(*pseudo_trials_t)))

    unit_code_inds = list(unit_grouper.groups.keys())
    unit_inds, unit_codes = list(map(list, list(zip(*[(ind[:len(units_index)], ind[-1])for ind in unit_code_inds]))))
    classes_unique = list(behavioral_units.groupby(class_condition_columns).groups.keys())
    classes_zips = ([(code_i, class_i)
                for code_i, class_i in enumerate(classes_unique)
                for _ in range(int(len(pseudo_trials) / len(classes_unique)))])
    class_codes, classes = list(map(list, zip(*classes_zips)))

    return pseudo_trials, unit_inds, unit_codes, class_codes, classes