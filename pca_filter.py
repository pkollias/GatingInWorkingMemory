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

    pca_condition_events_dict = md.np_loader(src_filename)

    unit_ind_list = []
    pca_matrix = []
    num_timebins = len(interval_split_to_bins_onset(t_start, t_end, timebin, timestep))
    unit_condition_fr_dims = (-1, num_timebins) ### TODO: correct for right number of rows? (len(condition_list), num_timebins)

    for unit_ind in pca_condition_events_dict.keys():
        if units.loc[unit_ind]['UnitNum'] != 0 and units.loc[unit_ind]['RatingCode'] != 7 and units.loc[unit_ind]['Area'] in area_list \
                and pca_condition_events_dict[unit_ind]['valid_conditions']:

            exceed_thr = lambda condition_fr_dict: len(condition_fr_dict) >= counts_thr
            if np.all(list(map(exceed_thr, pca_condition_events_dict[unit_ind]['unit_cond_dict'].values()))):

                unit_condition_fr = np.empty((0, num_timebins), float)

                for condition in condition_list:
                    condition_events_fr = list(pca_condition_events_dict[unit_ind]['unit_cond_dict'][condition].values())
                    counts = counts_thr if balance else len(condition_events_fr)
                    condition_events_counts_fr = sample(condition_events_fr, counts)
                    condition_fr = np.average(condition_events_counts_fr, axis=0).reshape((1, -1)) if average else np.array(condition_events_counts_fr)

                    unit_condition_fr = np.append(unit_condition_fr, condition_fr, axis=0)

                unit_condition_fr = unit_condition_fr.reshape((-1))
                unit_ind_list.append(unit_ind)
                pca_matrix.append(unit_condition_fr)

    pca_matrix = np.array(pca_matrix)
    conditions_events_filter_dict = {'pca_matrix': pca_matrix,
                                     'unit_condition_fr_dims': unit_condition_fr_dims,
                                     'unit_ind_list': unit_ind_list}

    md.np_saver(conditions_events_filter_dict, target_filename)


main()
