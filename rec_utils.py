from rec_analyses import *


def pbt_from_behavioral_units(condition_columns: list, version_fr: str, behavioral_units: pd.core.frame.DataFrame, db: DataBase) -> PopulationBehavioralTimeseries:

    v_fr_params = anova_version_fr_params(version_fr)
    t_start = v_fr_params['t_start']
    t_end = v_fr_params['t_end']
    timebin = v_fr_params['timebin']
    timestep = v_fr_params['timestep']
    # init
    timebin_interval = TimebinInterval(timebin, timestep, t_start, t_end)
    data_list = []
    # find units of behavioral units
    bu_grouper = behavioral_units.groupby(db.md.preproc_imports['units']['index'])
    # for every unit
    for unit_ind, unit_group in bu_grouper:
        # load firing rate data
        bufr = BehavioralUnitFiringRate(db, {'fr': version_fr}, unit_ind)
        # group by condition
        bu_cond_grouper = unit_group.groupby(condition_columns)
        # for every condition
        for condition, cond_group in bu_cond_grouper:
            # for every event
            for ii, (_, event_entry) in enumerate(cond_group.iterrows()):
                event_ind = tuple(event_entry[db.md.preproc_imports['events']['index']])
                condition = tuple(event_entry[condition_columns])
                timeseries = bufr.data[event_ind]
                data_entry = [unit_ind, event_ind, condition, ii, timeseries]
                # append to data_list
                data_list.append(data_entry)
    # create pbt
    pbt = PopulationBehavioralTimeseries(condition_columns, timebin_interval)
    pbt.add_data_rows_from_list(data_list)
    return pbt


def fbt_df_from_PCA(X_factor: np.ndarray, records, decomp_obj, timebin_interval):
    records_df = pd.DataFrame(records)
    num_factors = decomp_obj.n_components_
    num_conditions = records_df['Condition'].nunique()
    num_instances = records_df['Instance'].nunique()
    num_timebins = timebin_interval.num_of_bins()
    X_factor_base = X_factor.transpose().reshape((num_factors, num_conditions, num_instances, num_timebins))
    records_base = list(
        product(list(range(num_factors)), list(records_df['Condition'].unique()), list(range(num_instances))))
    records_base_df = pd.DataFrame(records_base, columns=['Factor', 'Condition', 'Instance'])
    records_base_df['Timeseries'] = list(X_factor_base.reshape(len(records_base_df), num_timebins))
    return records_base_df


def fbt_df_from_dPCA(X_factor: np.ndarray, records, decomp_obj, timebin_interval):
    records_df = pd.DataFrame(records)
    num_factors = decomp_obj.n_components
    num_conditions = records_df['Condition'].nunique()
    num_instances = records_df['Instance'].nunique()
    num_timebins = timebin_interval.num_of_bins()
    X_factor_base = X_factor.reshape((num_factors, num_conditions, num_instances, num_timebins))
    records_base = list(
        product(list(range(num_factors)), list(records_df['Condition'].unique()), list(range(num_instances))))
    records_base_df = pd.DataFrame(records_base, columns=['Factor', 'Condition', 'Instance'])
    records_base_df['Timeseries'] = list(X_factor_base.reshape(len(records_base_df), num_timebins))
    return records_base_df
