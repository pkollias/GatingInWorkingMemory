import sys
from rec import *
from pca_format import *
from random import sample


def main():


    # load analysis parameters
    args = sys.argv
    version_pca = args[1]
    version_fr = args[2]
    version_filter = args[3]
    version_units = args[4]
    counts_thr = int(args[5])



    # version parameters
    v_pca_params = pca_generate_conditions(version_pca)
    condition_list = v_pca_params['condition_list']

    v_fr_params = anova_version_fr_params(version_fr)
    t_start = v_fr_params['t_start']
    t_end = v_fr_params['t_end']
    timebin = v_fr_params['timebin']
    timestep = v_fr_params['timestep']

    v_filter_params = pca_version_filter_params(version_filter)
    balance = v_filter_params['balance']
    average = v_filter_params['average']
    if not (balance or average):
        print('Data needs to be either balanced or averaged', file=sys.stderr)
        exit()
    filter_params_list_str = v_filter_params['filter_params_list_str']

    v_units_area_list = pca_version_units_area_list(version_units)
    area_list = v_units_area_list['area_list']
    area_list_str = v_units_area_list['area_list_str']



    # load data and init vars, tables, and slices
    md = MetaData()
    db = md.db_base_loader(['units'])
    units = db['units']

    src_filename = md.proc_dest_path(path.join('BehavioralUnits', 'PCA', version_pca,
                                               behunit_params_str(version_fr, timebin, timestep, t_start, t_end), 'wrangle'),
                                     'conditions_events_dict.pkl')
    target_filename = md.proc_dest_path(path.join('BehavioralUnits', 'PCA', version_pca,
                                                  behunit_params_str(version_fr, timebin, timestep, t_start, t_end),
                                                  pca_filter_params_str(filter_params_list_str, counts_thr, area_list_str), 'filter'),
                                        'conditions_events_filter_dict.pkl')
    print(target_filename)
    if path.exists(target_filename):
        exit()

    # dictionary of unit_dicts -> dictionary of pca_conditions_dicts -> dictionary of events_lists/timeseries -> list of firing rates
    # aka timeseries for every event for every condition for every unit
    pca_condition_events_dict = md.np_loader(src_filename)

    unit_ind_list = []
    pca_matrix = []
    num_timebins = len(interval_split_to_bins_onset(t_start, t_end, timebin, timestep))
    unit_condition_fr_dims = (-1, num_timebins) ### TODO: correct for right number of rows? (len(condition_list), num_timebins)

    # for every unit
    for unit_ind in pca_condition_events_dict.keys():

        # if not multiunit, belongs in area of interest, and has at least one event from all pca_conditions
        if units.loc[unit_ind]['UnitNum'] != 0 and units.loc[unit_ind]['RatingCode'] != 7 and units.loc[unit_ind]['Area'] in area_list \
                and pca_condition_events_dict[unit_ind]['valid_conditions']:

            exceed_thr = lambda condition_fr_dict: len(condition_fr_dict) >= counts_thr
            # if all pca_conditions have minimum number of events
            if np.all(list(map(exceed_thr, pca_condition_events_dict[unit_ind]['unit_cond_dict'].values()))):

                # firing rate template for timeseries
                unit_condition_fr = np.empty((0, num_timebins), float)

                # for every pca_condition
                for condition in condition_list:
                    # get list of all timeseries for that pca_condition for that unit
                    condition_events_fr = list(pca_condition_events_dict[unit_ind]['unit_cond_dict'][condition].values())
                    # set number of timeseries that will go into finalized table for that pca_condition (threshold if balanced or max otherwise)
                    counts = counts_thr if balance else len(condition_events_fr)
                    # counts x num_timebins sized list of timeseries
                    condition_events_counts_fr = sample(condition_events_fr, counts)
                    # convert to array and average if necessary
                    condition_fr = np.average(condition_events_counts_fr, axis=0).reshape((1, -1)) if average else np.array(condition_events_counts_fr)
                    # append to list of pca_conditions
                    # if not averaged list will be (num_pca_conditions * counts) x num_timebins with each pca_condition blocked in sets of counts
                    unit_condition_fr = np.append(unit_condition_fr, condition_fr, axis=0)

                # reshape to (num_pca_conditions * counts * num_timebins) x 1
                # where each condition is blocked in sets of instances where each instance is blokced in its timeseries
                # [pca_cond_1_instance_1_t0, pca_cond_1_instance_1_t1, ..., pca_cond_1_instance_2_t0, ..., pca_cond_2_instance_1_t0, ...]
                unit_condition_fr = unit_condition_fr.reshape((-1))
                # append to list of units
                unit_ind_list.append(unit_ind)
                pca_matrix.append(unit_condition_fr)

    pca_matrix = np.array(pca_matrix)
    conditions_events_filter_dict = {'pca_matrix': pca_matrix,
                                     'unit_ind_list': unit_ind_list}

    md.np_saver(conditions_events_filter_dict, target_filename)


main()
