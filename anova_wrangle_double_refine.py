import sys
from rec import *
from functools import reduce
from rec_format import *
from itertools import product
from versioning import *


def main():

    # load analysis parameters
    args = sys.argv
    u_iloc = int(args[1])
    version_aov_list = args[2]
    version_fr = args[3]

    valid_thr = 15

    # version parameters
    anova_selection_dict_dict = {}
    src_filename_dict = {}
    # load all anova_dicts
    for version_aov in version_aov_list.split('_'):

        v_fr_params = anova_version_fr_params(version_fr)
        t_start = v_fr_params['t_start']
        t_end = v_fr_params['t_end']
        timebin = v_fr_params['timebin']
        timestep = v_fr_params['timestep']

        # load data and init vars, tables, and slices
        md = MetaData()
        db = md.db_base_loader(['trials', 'events', 'units', 'conditions'])
        trials, events, units, conditions = db['trials'], db['events'], db['units'], db['conditions']
        unit_entry = units.iloc[u_iloc]
        unit_ind = tuple(unit_entry[md.preproc_imports['units']['index']])
        sess, channum, unitnum = unit_ind
        src_filename_dict[version_aov] = md.proc_dest_path(path.join('BehavioralUnits', 'Anova', version_aov,
                                                                     behunit_params_str(version_fr, timebin, timestep, t_start, t_end), 'wrangle'),
                                                           'selection_dict_{0:s}_chan{1:03d}_unit{2:03d}.pkl'.
                                                           format(sess, channum, unitnum))

        anova_selection_dict_dict[version_aov] = md.np_loader(src_filename_dict[version_aov])


    # evaluate anova_dicts
    valid_list = [anova_selection_dict['valid'] for anova_selection_dict in anova_selection_dict_dict.values()]
    num_events_list = [anova_selection_dict['num_events'] for anova_selection_dict in anova_selection_dict_dict.values()]
    keys_list = list(anova_selection_dict_dict.keys())
    if not all(valid_list):
        for anova_selection_dict in anova_selection_dict_dict.values():
            anova_selection_dict['df'] = np.nan
            anova_selection_dict['selection_df'] = np.nan
            anova_selection_dict['counts'] = np.nan
            anova_selection_dict['valid'] = False
            anova_selection_dict['num_events'] = 0
    else:
        min_val = min(num_events_list)
        max_index = np.argmax(num_events_list)

        version_aov = keys_list[max_index]
        anova_selection_dict = anova_selection_dict_dict[version_aov]
        v_aov_params = anova_version_aov_params(version_aov, version_fr)
        selection_dict = v_aov_params['selection_dict']
        levels_dict = v_aov_params['levels_dict']

        # split dfs into selection dfs and subsample min number events
        level_combinations = list(product(*list(levels_dict.values())))
        df = anova_selection_dict['df']
        counts = anova_selection_dict['counts']
        selection_df = {}
        for selection in selection_dict['list']:
            full_selection_df = df.loc[selection].reset_index(drop=True)
            level_comb_series = full_selection_df.apply(lambda row: tuple(row[x] for x in levels_dict.keys()), axis=1)
            selection_df[selection] = pd.concat([full_selection_df.loc[level_comb_series == level].sample(n=min_val, random_state=0)
                                                 for level in level_combinations], axis=0)

        # anova selection_dict create and overwrite
        anova_selection_dict = {'df': df,
                                'selection_df': selection_df,
                                'counts': counts,
                                'num_events': min_val,
                                'valid': True}
        md.np_saver(anova_selection_dict, src_filename_dict[version_aov])


main()
