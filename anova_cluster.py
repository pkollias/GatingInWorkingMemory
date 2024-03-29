import sys
from rec import *
from aov_stats import *
from versioning import *


def main():
    """ u_iloc, aov, fr, selection, shuffles, [overwrite] """
    args_version = sys.argv[1:]
    # args_version = ['job_id=0', 'overwrite=True']
    version = job_scheduler(args_version, args_from_parse_func)

    u_iloc = int(version['u_iloc'])
    version_aov = version['aov']
    selection = version['selection']
    shuffles = int(version['shuffles'])
    version_fr = version['fr']

    omega_percentile = 95
    cluster_percentile = 95

    # version parameters
    v_aov_params = anova_version_aov_params(version_aov, version_fr)
    x_factors = v_aov_params['x_factors']

    v_fr_params = version_fr_params(version_fr)
    t_start = v_fr_params['t_start']
    t_end = v_fr_params['t_end']
    timebin = v_fr_params['timebin']
    timestep = v_fr_params['timestep']



    # load data and init vars, tables, and slices
    md = MetaData()
    db = md.db_base_loader(['units'])
    units = db['units']
    unit_entry = units.iloc[u_iloc]
    sess, channum, unitnum = unit_entry.name
    src_filename = md.proc_dest_path(path.join('BehavioralUnits', 'Anova', version_aov,
                                               behunit_params_str(version_fr, timebin, timestep, t_start, t_end),
                                               'time_anova', selection),
                                     'time_anova_results_{0:s}_chan{1:03d}_unit{2:03d}.pkl'.
                                     format(sess, channum, unitnum))
    target_filename = md.proc_dest_path(path.join('BehavioralUnits', 'Anova', version_aov,
                                                  behunit_params_str(version_fr, timebin, timestep, t_start, t_end),
                                                  'time_cluster', selection),
                                        'time_anova_cluster_results_{0:s}_chan{1:03d}_unit{2:03d}.pkl'.
                                        format(sess, channum, unitnum))
    print(target_filename)
    if path.exists(target_filename) and ('overwrite' not in version.keys() or not eval(version['overwrite'])):
        exit()

    # load unit time results
    unit_time_results = md.np_loader(src_filename)

    unit_time_cluster_results = {}
    # for every factor in anova
    for x_i in x_factors:

        # initialize cluster results
        unit_time_cluster_results[x_i] = {}
        # round anova results, correct all 0 anovas and return omega sq difference from threshold, threshold and zscore
        omega_dict = dict([(int(y_bin_str_to_strnum(bin)),
                            evaluate_anova_shuffles(results['results']['shuffles'][x_i], omega_percentile, shuffles))
                           for bin, results in unit_time_results.items()])

        # initialize distribution of cluster values
        cluster_val_distr = []
        # for every shuffle
        for shuffle_i in range(1, shuffles):
            # get list of omega differences from anova result dict
            omega_sq_diff_list = [(key, val['omega_sq_diff'].loc[shuffle_to_name(shuffle_i)]) for key, val in omega_dict.items()]
            # evaluate differences and get list of clusters with values and timebins
            shuffle_clusters = get_clusters(clusters_from_shuffle_list(omega_sq_diff_list))
            if shuffle_clusters:
                # append maximum value cluster
                cluster_val_distr.append(max([cluster['val'] for cluster in shuffle_clusters]))
            else: ### TODO: Is that the right way to do cluster correction?
                cluster_val_distr.append(0)

        # estimate cluster threshold
        cluster_threshold = np.percentile(cluster_val_distr, cluster_percentile, interpolation='linear') if cluster_val_distr else -np.inf

        # repeat process for observed clusters
        observed_omega_sq_diff_list = [(key, val['omega_sq_diff'].loc[shuffle_to_name(0)]) for key, val in omega_dict.items()]
        observed_clusters = get_clusters(clusters_from_shuffle_list(observed_omega_sq_diff_list))
        # evaluate observed significant clusters
        significant_shuffle_clusters = [cluster['bins'] for cluster in observed_clusters if cluster['val'] > cluster_threshold]

        unit_time_cluster_results[x_i] = {'omega': omega_dict,
                                          'cluster_val_distr': cluster_val_distr,
                                          'cluster_threshold': cluster_threshold,
                                          'clusters': significant_shuffle_clusters}

    md.np_saver(unit_time_cluster_results, target_filename)


def args_from_parse_func(parse_version):

    args_version_list = []

    args_u_iloc = ['u_iloc={0:d}'.format(u_iloc) for u_iloc in range(2436)]
    args_aov = ['aov={0:s}'.format(aov) for aov in ['PresentedStimulus', 'GatedStimulus']]
    args_selection = ['selection={0:s}'.format(selection) for selection in ['Cue', 'PreDist', 'Gating', 'PostDist', 'Target']]
    args_shuffles = ['shuffles=2000']
    args_fr = ['fr=ConcatFactor']
    args_version_list.extend(list(map(list, list(product(args_u_iloc, args_aov, args_selection, args_shuffles, args_fr)))))

    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))

    return args_version_from_job


main()
