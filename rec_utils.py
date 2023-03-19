from rec_analyses import *
from functools import reduce


def pbt_from_behavioral_units(condition_columns: list, version_fr: str, behavioral_units: pd.core.frame.DataFrame, db: DataBase) -> PopulationBehavioralTimeseries:

    v_fr_params = version_fr_params(version_fr)
    t_start = v_fr_params['t_start']
    t_end = v_fr_params['t_end']
    timebin = v_fr_params['timebin']
    timestep = v_fr_params['timestep']
    md = db.md
    units_index, events_index = md.preproc_imports['units']['index'], md.preproc_imports['events']['index']
    # init
    timebin_interval = TimebinInterval(timebin, timestep, t_start, t_end)
    data_list = []
    # find units of behavioral units
    bu_grouper = behavioral_units.groupby(units_index + ['Unit_Code'])
    # for every unit
    for unit_ind_code, unit_group in bu_grouper:
        # load firing rate data
        unit_ind = unit_ind_code[:len(units_index)]
        unit_code = unit_ind_code[-1]
        bufr = BehavioralUnitFiringRate(db, {'fr': version_fr}, unit_ind)
        # group by condition
        bu_cond_grouper = unit_group.groupby(condition_columns)
        print(unit_ind_code)
        # for every condition
        for condition, cond_group in bu_cond_grouper:
            # for every event
            for ii, (_, event_entry) in enumerate(cond_group.iterrows()):
                event_ind = tuple(event_entry[db.md.preproc_imports['events']['index']])
                condition = tuple(event_entry[condition_columns])
                timeseries = bufr.data[event_ind]
                data_entry = [unit_ind, unit_code, event_ind, condition, ii, timeseries]
                # append to data_list
                data_list.append(data_entry)
    # create pbt
    pbt = PseudoPopulationBehavioralTimeseries(condition_columns, timebin_interval)
    pbt.add_data_rows_from_list(data_list)
    return pbt


def fbt_df_from_PCA(X_factor: np.ndarray, records, num_factors, timebin_interval):
    records_df = pd.DataFrame(records)
    num_conditions = records_df['Condition'].nunique()
    num_instances = records_df['Instance'].nunique()
    num_timebins = timebin_interval.num_of_bins()
    X_factor_base = X_factor.transpose().reshape((num_factors, num_conditions, num_instances, num_timebins))
    records_base = list(
        product(list(range(num_factors)), list(records_df['Condition'].unique()), list(range(num_instances))))
    records_base_df = pd.DataFrame(records_base, columns=['Factor', 'Condition', 'Instance'])
    records_base_df['Timeseries'] = list(X_factor_base.reshape(len(records_base_df), num_timebins))
    return records_base_df


def fbt_df_from_dPCA(X_factor: np.ndarray, records, num_factors, timebin_interval):
    records_df = pd.DataFrame(records)
    num_conditions = records_df['Condition'].nunique()
    num_instances = records_df['Instance'].nunique()
    num_timebins = timebin_interval.num_of_bins()
    X_factor_base = X_factor.reshape((num_factors, num_conditions, num_instances, num_timebins))
    records_base = list(
        product(list(range(num_factors)), list(records_df['Condition'].unique()), list(range(num_instances))))
    records_base_df = pd.DataFrame(records_base, columns=['Factor', 'Condition', 'Instance'])
    records_base_df['Timeseries'] = list(X_factor_base.reshape(len(records_base_df), num_timebins))
    return records_base_df


def conjoint_units_behavioral_lists(units_behavioral_lists_list, filename_list):

    md = MetaData()

    index = units_behavioral_lists_list[0].index
    mask = reduce(lambda x, y: ~(x.apply(bool) & y.apply(bool)), units_behavioral_lists_list)
    empty = pd.Series([[] for _ in range(len(index))], index=index)

    for valid_units_behavioral_lists, target_filename in zip(units_behavioral_lists_list, filename_list):
        md.np_saver(valid_units_behavioral_lists.mask(mask, empty), target_filename)


def fr_from_units_events(units_events_filter: pd.core.frame.DataFrame, version_fr: str, db: DataBase, t_ind: int) -> np.ndarray:

    units_events_fr = units_events_filter.drop('valid', axis=1)

    for unit_ind in units_events_fr.columns:

        bufr = BehavioralUnitFiringRate(db, {'fr': version_fr}, unit_ind)
        units_events_fr[unit_ind] = pd.Series(bufr.data).apply(lambda x: x[t_ind]).reindex(units_events_fr.index)

    return units_events_fr
