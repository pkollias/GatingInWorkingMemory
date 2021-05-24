import sys
from rec import *
from aov_stats import *
from versioning import *


def main():

    # load analysis parameters
    args = sys.argv
    version_aov = args[1]
    version_fr = args[2]
    shuffles = int(args[3])


    omega_percentile = 95

    # version parameters
    v_aov_params = anova_version_aov_params(version_aov, version_fr)
    x_factors = v_aov_params['x_factors']
    selection_dict = v_aov_params['selection_dict']

    v_fr_params = anova_version_fr_params(version_fr)
    t_start = v_fr_params['t_start']
    t_end = v_fr_params['t_end']
    timebin = v_fr_params['timebin']
    timestep = v_fr_params['timestep']



    # load data and init vars, tables, and slices
    md = MetaData()
    db = md.db_base_loader(['physiology'])
    physiology = db['physiology']
    target_filename = md.proc_dest_path(path.join('BehavioralUnits', 'Anova', version_aov,
                                                  behunit_params_str(version_fr, timebin, timestep, t_start, t_end), 'summarize'),
                                        'physiology_dict.pkl')
    print(target_filename)
    if path.exists(target_filename):
        exit()

    physiology_dict = {}
    # for every selection
    for selection in selection_dict['list']:

        # init selection dict
        physiology_dict[selection] = {}
        # for every unit
        for ind, unit_index in enumerate(physiology.index):

            print(selection, ind, unit_index)
            sess, channum, unitnum = unit_index
            src_filename = md.proc_dest_path(path.join('BehavioralUnits', 'Anova', version_aov,
                                                       behunit_params_str(version_fr, timebin, timestep, t_start,
                                                                          t_end),
                                                       'time_anova', selection),
                                             'time_anova_results_{0:s}_chan{1:03d}_unit{2:03d}.pkl'.
                                             format(sess, channum, unitnum))
            # load anova unit time results
            unit_time_results = md.np_loader(src_filename)
            # init unit dict
            physiology_dict[selection][unit_index] = {}

            omega_dict = {}
            for x_i in x_factors:

                # init x_factor dict
                valid = bool(unit_time_results)
                significant = []
                zscores = []

                if valid:

                    # round anova results, correct all 0 anovas and return omega sq difference from threshold, threshold and zscore
                    omega_dict = dict([(int(y_bin_str_to_strnum(bin)),
                                        evaluate_anova_shuffles(results['results']['shuffles'][x_i], omega_percentile,
                                                                shuffles))
                                       for bin, results in unit_time_results.items()])
                    # add significant zscores and significant timebins
                    significant = [key for key, val in omega_dict.items() if val['omega_sq_observed'] > val['threshold']]
                    zscores = [val['zscore'] for val in omega_dict.values()]

                physiology_dict[selection][unit_index][x_i] = {'valid': valid,
                                                               'significant': significant,
                                                               'zscores': zscores}


    md.np_saver(physiology_dict, target_filename)


main()
