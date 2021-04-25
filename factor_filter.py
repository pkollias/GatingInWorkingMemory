import sys
from rec import *
from rec_stats import *
from rec_format import *
from random import sample, seed


def main():


    # load analysis parameters
    args = sys.argv
    version_factor = args[1]
    version_fr = args[2]
    version_filter = args[3]
    version_areas = args[4]
    version_subjects = args[5]
    counts_thr = int(args[6])
    seed_str = '_'.join(args[1:6])


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
    filter_params_list_str = v_filter_params['filter_params_list_str']

    v_units_area_list = factor_version_units_area_list(version_areas)
    area_list = v_units_area_list['area_list']
    area_list_str = v_units_area_list['area_list_str']

    v_units_subject_list = factor_version_units_subject_list(version_subjects)
    subject_list = v_units_subject_list['subject_list']
    subject_list_str = v_units_subject_list['subject_list_str']

    units_str = '_'.join([area_list_str, subject_list_str])

    exceed_thr = lambda condition_fr_dict: len(condition_fr_dict) >= counts_thr
    seed(seed_str)

    event_fr_threshold = 100



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

    timebin_interval = TimebinInterval(timebin, timestep, t_start, t_end)
    pbt = PopulationBehavioralTimeseries(condition_columns, timebin_interval)
    data_list = []
    # for every unit
    for unit_ind in factor_condition_events_dict.keys():

        # if not multiunit, belongs in area of interest, at least one event from all factor_conditions and has min count of events
        if units.loc[unit_ind]['UnitNum'] != 0 and units.loc[unit_ind]['RatingCode'] != 7 and units.loc[unit_ind]['Area'] in area_list \
                and factor_condition_events_dict[unit_ind]['valid_conditions'] \
                and sessions.loc[units.loc[unit_ind].Session].Subject in subject_list \
                and np.all(list(map(exceed_thr, factor_condition_events_dict[unit_ind]['unit_cond_dict'].values()))):

            # load unit's firing rate data
            src_filename = md.proc_dest_path(path.join('BehavioralUnits', 'FiringRates',
                                                       behunit_params_str(version_fr, timebin, timestep, t_start,
                                                                          t_end)),
                                             'behunit_FR_{0:s}_chan{1:03d}_unit{2:03d}.pkl'.format(*unit_ind))
            timebin_fr_dict = md.np_loader(src_filename)

            # threshold for events with abnormally high firing rate
            if np.all([exceed_thr(list(filter(lambda event_ind: np.mean(timebin_fr_dict[event_ind]) < event_fr_threshold, event_index_list)))
                       for event_index_list
                       in factor_condition_events_dict[unit_ind]['unit_cond_dict'].values()]):

                # for every factor_condition
                for condition in condition_list:
                    # get list of all event indices for that factor_condition for that unit
                    condition_events_tuples = [(event_ind, timebin_fr_dict[event_ind])
                                               for event_ind
                                               in factor_condition_events_dict[unit_ind]['unit_cond_dict'][condition]
                                               if np.mean(timebin_fr_dict[event_ind]) < event_fr_threshold]
                    # set number of timeseries that will go into finalized table for that factor_condition
                    counts = counts_thr
                    # (1 x counts of inds, counts x num_timebins of timeseries)
                    condition_events_counts_tuples = sample(condition_events_tuples, counts)
                    for inst_ind, (event_ind, timeseries) in enumerate(condition_events_counts_tuples):
                        data_list.append([unit_ind, event_ind, condition, inst_ind, timeseries])

    pbt.add_data_rows_from_list(data_list)


    md.np_saver(pbt, target_filename)


main()
