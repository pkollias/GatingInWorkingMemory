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

    # version parameters
    v_aov_params = anova_version_aov_params(version_aov, version_fr)
    levels_dict = v_aov_params['levels_dict']
    group_column_list = v_aov_params['group_column_list']

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
                                               'wrangle'),
                                     'selection_dict_{0:s}_chan{1:03d}_unit{2:03d}.pkl'.
                                     format(sess, channum, unitnum))
    target_filename = md.proc_dest_path(path.join('BehavioralUnits', 'Anova', version_aov,
                                                  behunit_params_str(version_fr, timebin, timestep, t_start, t_end),
                                                  'time_anova', selection),
                                        'time_anova_results_{0:s}_chan{1:03d}_unit{2:03d}.pkl'.
                                        format(sess, channum, unitnum))
    print(target_filename)
    if path.exists(target_filename) and ('overwrite' not in version.keys() or not eval(version['overwrite'])):
        exit()

    events_index = md.preproc_imports['events']['index']
    x_list = list(levels_dict.keys())

    # load anova selection dict
    anova_selection_dict = md.np_loader(src_filename)

    unit_time_results = {}
    # if valid unit
    if anova_selection_dict['valid']:

        # select df with selection subset of events
        df = anova_selection_dict['selection_df'][selection]
        # select timebin columns
        y_strnum_list = list(df.drop(set(events_index + x_list + group_column_list), axis=1).columns) # TODO: change wrangle column names
        # convert timebin column names into formatted strings
        y_str_list = [y_bin_strnum_to_str(y_strnum) for y_strnum in y_strnum_list]
        y_rename_dict =dict(zip(y_strnum_list, y_str_list))
        df.rename(columns=y_rename_dict, inplace=True)
        # for every timebin
        for y_str, y_numstr in zip(y_str_list, y_strnum_list):
            # if all values the same only run observed anova
            num_shuffles = 1 if len(df[y_str].unique()) == 1 else shuffles
            # anova results
            unit_results = aov_shuffle_and_results(df, y_str, x_list, num_shuffles, group_column_list)
            unit_time_results[y_str] = {'results': unit_results,
                                        'bin_onset': int(y_numstr)}

    md.np_saver(unit_time_results, target_filename)


def args_from_parse_func(parse_version):

    args_version_list = []

    args_u_iloc = ['u_iloc={0:d}'.format(u_iloc) for u_iloc in range(2436)]
    args_aov = ['aov={0:s}'.format(aov) for aov in ['GatedCue']]
    args_selection = ['selection={0:s}'.format(selection) for selection in ['Cue', 'PreDist', 'Gating', 'PostDist', 'Target']]
    args_shuffles = ['shuffles=2000']
    args_fr = ['fr=ConcatFactor']
    args_version_list.extend(list(map(list, list(product(args_aov, args_selection, args_shuffles, args_fr, args_u_iloc)))))

    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))

    return args_version_from_job


main()
