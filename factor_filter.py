import sys
from rec import *
from rec_stats import *
from factor_format import *
from random import sample


def main():


    # load analysis parameters
    args = sys.argv
    version_factor = args[1]
    version_fr = args[2]
    version_filter = args[3]
    version_areas = args[4]
    version_subjects = args[5]
    counts_thr = int(args[6])



    # version parameters
    v_factor_params = factor_generate_conditions(version_factor)
    condition_columns = v_factor_params['condition_columns']
    condition_list = v_factor_params['condition_list']

    v_fr_params = anova_version_fr_params(version_fr)
    t_start = v_fr_params['t_start']
    t_end = v_fr_params['t_end']
    timebin = v_fr_params['timebin']
    timestep = v_fr_params['timestep']

    v_filter_params = factor_version_filter_params(version_filter)
    balance = v_filter_params['balance']
    average = v_filter_params['average']
    if not (balance or average):
        print('Data needs to be either balanced or averaged', file=sys.stderr)
        exit()
    filter_params_list_str = v_filter_params['filter_params_list_str']

    v_units_area_list = factor_version_units_area_list(version_areas)
    area_list = v_units_area_list['area_list']
    area_list_str = v_units_area_list['area_list_str']

    v_units_subject_list = factor_version_units_subject_list(version_subjects)
    subject_list = v_units_subject_list['subject_list']
    subject_list_str = v_units_subject_list['subject_list_str']

    units_str = '_'.join([area_list_str, subject_list_str])



    # load data and init vars, tables, and slices
    md = MetaData()
    db = md.db_base_loader(['units', 'sessions'])
    units, sessions = db['units'], db['sessions']

    src_filename = md.proc_dest_path(path.join('BehavioralUnits', 'Factorization', version_factor,
                                               behunit_params_str(version_fr, timebin, timestep, t_start, t_end), 'wrangle'),
                                     'conditions_events_dict.pkl')
    target_filename = md.proc_dest_path(path.join('BehavioralUnits', 'Factorization', version_factor,
                                                  behunit_params_str(version_fr, timebin, timestep, t_start, t_end),
                                                  factor_filter_params_str(filter_params_list_str, counts_thr, units_str), 'filter'),
                                        'conditions_events_filter_obj.pkl')
    print(target_filename)
    if path.exists(target_filename):
        exit()

    # dictionary of unit_dicts -> dictionary of factor_conditions_dicts -> dictionary of events_lists/timeseries -> list of firing rates
    # aka timeseries for every event for every condition for every unit
    factor_condition_events_dict = md.np_loader(src_filename)

    num_timebins = len(interval_split_to_bins_onset(t_start, t_end, timebin, timestep))
    num_instances = 1 if average else counts_thr
    num_conditions = len(condition_list)
    unit_condition_fr_dims = (-1, num_timebins)
    # unit_ind_list = []
    # factor_matrix = []

    pbt = PopulationBehavioralTimeseries(shape=(0, num_conditions, num_instances, num_timebins), condition_labels=condition_columns,
                                         timebins={'version_fr': version_fr, 'timebin': timebin, 'timestep': timestep})

    # for every unit
    for unit_ind in factor_condition_events_dict.keys():

        # if not multiunit, belongs in area of interest, and has at least one event from all factor_conditions
        if units.loc[unit_ind]['UnitNum'] != 0 and units.loc[unit_ind]['RatingCode'] != 7 and units.loc[unit_ind]['Area'] in area_list \
                and factor_condition_events_dict[unit_ind]['valid_conditions'] \
                and sessions.loc[units.loc[unit_ind].Session].Subject in subject_list:

            exceed_thr = lambda condition_fr_dict: len(condition_fr_dict) >= counts_thr
            # if all factor_conditions have minimum number of events
            if np.all(list(map(exceed_thr, factor_condition_events_dict[unit_ind]['unit_cond_dict'].values()))):

                # create behavioral timeseries instance for current unit
                ubt = pbt.derive_unit(unit_ind)

                # for every factor_condition
                for condition in condition_list:
                    # get list of all timeseries for that factor_condition for that unit
                    condition_events_fr = list(factor_condition_events_dict[unit_ind]['unit_cond_dict'][condition].values())
                    # set number of timeseries that will go into finalized table for that factor_condition (threshold if balanced or max otherwise)
                    counts = counts_thr if balance else len(condition_events_fr)
                    # counts x num_timebins sized list of timeseries
                    condition_events_counts_fr = sample(condition_events_fr, counts)
                    # convert to array and average if necessary
                    condition_fr = np.average(condition_events_counts_fr, axis=0).reshape((1, -1)) if average else np.array(condition_events_counts_fr)

                    # uct_shape = (counts, num_timebins)
                    # create condition timeseries instance for current unit current condition
                    uct = ubt.derive_unit_condition(condition_levels=condition)
                    uct.set_base_data(condition_fr)

                    ubt.add_condition(uct)

                pbt.add_unit(ubt)



    md.np_saver(pbt, target_filename)


main()
