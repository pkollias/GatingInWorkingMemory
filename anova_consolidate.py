import sys
from rec import *
from rec_format import *
from versioning import *


def main():
    """ aov, fr, [overwrite] """
    args_version = sys.argv[1:]
    # args_version = ['job_id=0', 'overwrite=True']
    version = job_scheduler(args_version, args_from_parse_func)

    version_aov = version['aov']
    version_fr = version['fr']

    # version parameters
    v_aov_params = anova_version_aov_params(version_aov, version_fr)
    selection_dict = v_aov_params['selection_dict']

    v_fr_params = version_fr_params(version_fr)
    t_start = v_fr_params['t_start']
    t_end = v_fr_params['t_end']
    timebin = v_fr_params['timebin']
    timestep = v_fr_params['timestep']



    # load data and init vars, tables, and slices
    md = MetaData()
    db = md.db_base_loader(['physiology'])
    physiology = db['physiology']
    target_filename = md.proc_dest_path(path.join('BehavioralUnits', 'Anova', version_aov,
                                                  behunit_params_str(version_fr, timebin, timestep, t_start, t_end), 'consolidate'),
                                        'physiology_dict.pkl')
    print(target_filename)
    if path.exists(target_filename) and ('overwrite' not in version.keys() or not eval(version['overwrite'])):
        exit()

    physiology_dict = {}

    for selection in selection_dict['list']:

        physiology_dict[selection] = {}
        for ind, unit_index in enumerate(physiology.index):

            print(selection, ind, unit_index)

            sess, channum, unitnum = unit_index
            src_filename = md.proc_dest_path(path.join('BehavioralUnits', 'Anova', version_aov,
                                                       behunit_params_str(version_fr, timebin, timestep, t_start, t_end),
                                                       'time_cluster', selection),
                                             'time_anova_cluster_results_{0:s}_chan{1:03d}_unit{2:03d}.pkl'.
                                             format(sess, channum, unitnum))
            unit_time_cluster_results = md.np_loader(src_filename)

            physiology_dict[selection][unit_index] = {}
            for x_i, cluster_results in unit_time_cluster_results.items():

                valid = bool(unit_time_cluster_results[x_i]['omega'])
                clusters = np.nan
                zscores = np.nan
                if valid:
                    clusters = unit_time_cluster_results[x_i]['clusters']
                    zscores = [omega['zscore'] for omega in unit_time_cluster_results[x_i]['omega'].values()]

                physiology_dict[selection][unit_index][x_i] = {'valid': valid,
                                                               'clusters': clusters,
                                                               'zscores': zscores}

    md.np_saver(physiology_dict, target_filename)


def args_from_parse_func(parse_version):

    args_version_list = []

    args_aov = ['aov={0:s}'.format(aov) for aov in ['GatedCue']]
    args_fr = ['fr=ConcatFactor']
    args_version_list.extend(list(map(list, list(product(args_aov, args_fr)))))

    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))

    return args_version_from_job


main()
