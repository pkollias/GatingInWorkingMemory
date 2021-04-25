import sys
from rec import *
from functools import reduce
from rec_format import *
from itertools import product


def main():

    ### TODO: class for firing rate versioning (firing rates for units, interval steps, times, ...)
    ### TODO: class for dfs and filtering
    ### TODO: class for all steps of anova (analysis version, selections, stats, intermediate results, clusters, ...)

    # load analysis parameters
    args = sys.argv
    u_iloc = int(args[1])
    version_aov = args[2]
    version_fr = args[3]


    valid_thr = 15

    # version parameters
    v_aov_params = anova_version_aov_params(version_aov, version_fr)
    selection_dict = v_aov_params['selection_dict']
    levels_dict = v_aov_params['levels_dict']
    event_cnj_mask = v_aov_params['event_cnj_mask']
    group_column_list = v_aov_params['group_column_list']

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
    src_filename = md.proc_dest_path(path.join('BehavioralUnits', 'FiringRates',
                                               behunit_params_str(version_fr, timebin, timestep, t_start, t_end)),
                                     'behunit_FR_{0:s}_chan{1:03d}_unit{2:03d}.pkl'.format(sess, channum, unitnum))
    target_filename = md.proc_dest_path(path.join('BehavioralUnits', 'Anova', version_aov,
                                                  behunit_params_str(version_fr, timebin, timestep, t_start, t_end), 'wrangle'),
                                        'selection_dict_{0:s}_chan{1:03d}_unit{2:03d}.pkl'.
                                        format(sess, channum, unitnum))
    print(target_filename)
    if path.exists(target_filename):
        exit()

    # load unit firing rate data
    timebin_fr_dict = md.np_loader(src_filename)

    # Create timebin_fr_df
    sm = SamplingMethod()
    event_window = SamplingInterval(SamplingPoint(sm, t=t_start * qu.ms),
                                    SamplingPoint(sm, t=t_end * qu.ms))
    event_fr_time_bins = event_window.bins_by_step(SamplingPoint(sm, stamp=sm.stamp_from_t(timebin * qu.ms)),
                                                   SamplingPoint(sm, t=timestep * qu.ms))
    bin_onset_list = [int(sm.t_from_stamp(event_fr_time_bin.start.stamp).to('ms').value) for event_fr_time_bin in event_fr_time_bins]
    bin_onset_str_list = [str(bin_onset) for bin_onset in bin_onset_list]
    timebin_fr_df = pd.DataFrame.from_dict(timebin_fr_dict, orient='index', columns=bin_onset_str_list)

    trials_index = md.preproc_imports['trials']['index']
    events_index = md.preproc_imports['events']['index']
    columns_trials = trials_index
    columns_valid = [selection_dict['column']] + list(levels_dict.keys())
    columns_anova = columns_valid + group_column_list
    columns_events_conditions = events_index + columns_anova
    mask_and = lambda mask_i, mask_j: mask_i & mask_j

    # filter correct only trials
    trials_slice = trials[trials['StopCondition'].eq(1) & trials['Session'].eq(sess)].reset_index(drop=True)[columns_trials]

    # filter anova events data (set of all selections and list of all levels from all factors)
    events_slice = events.loc[sess].reset_index(drop=True)
    conditions_slice = conditions.loc[conditions['Session'].eq(sess)]
    events_conditions = pd.merge(events_slice, conditions_slice, on=events_index)
    # calculate masks
    mask_true = pd.Series(True, index=np.arange(len(events_conditions)))
    selection_mask = mask_true if not bool(selection_dict) else events_conditions[selection_dict['column']].isin(selection_dict['list'])
    levels_mask = reduce(mask_and, [events_conditions[x].isin(levels) for x, levels in levels_dict.items()], mask_true)
    cnj_mask = reduce(mask_and, [filter_df_wrapper(events_conditions, mask_dict['column'], mask_dict['wrapper'], mask_dict['arg'])['mask']
                                 for mask_dict in event_cnj_mask], mask_true)
    ec_mask = (selection_mask & levels_mask & cnj_mask)
    events_conditions_slice = events_conditions.loc[ec_mask][columns_events_conditions]

    # merge into final events_conditions_df
    events_conditions_df = pd.merge(events_conditions_slice, trials_slice, on=trials_index)
    for x_i, levels_i in levels_dict.items():
        events_conditions_df[x_i].cat.set_categories(levels_i, inplace=True)
    events_conditions_df.set_index(events_index, inplace=True, drop=False)

    # create final df by concatanating events_conditions_df and timebin_fr_df
    df = pd.concat([events_conditions_df, timebin_fr_df], join='inner', axis=1).sort_index()
    # counts of events table
    counts = df.groupby(columns_valid).size()
    # minimum number of events in table
    selection_level_combinations = list(product(*[selection_dict['list']] + list(levels_dict.values())))
    all_combinations_in_df = not (bool(set(counts.index).difference(set(selection_level_combinations))) or
                                  bool(set(selection_level_combinations).difference(set(counts.index))))
    num_events = df.groupby(columns_valid).size().min()
    # if valid unit (all level combination values and exceeding min threshold)
    valid = all_combinations_in_df and num_events >= valid_thr

    # split dfs into selection dfs and subsample min number events
    level_combinations = list(product(*list(levels_dict.values())))
    df.set_index(selection_dict['column'], inplace=True)

    selection_df = {}
    if valid:
        for selection in selection_dict['list']:
            full_selection_df = df.loc[selection].reset_index(drop=True)
            level_comb_series = full_selection_df.apply(lambda row: tuple(row[x] for x in levels_dict.keys()), axis=1)
            selection_df[selection] = pd.concat([full_selection_df.loc[level_comb_series == level].sample(n=num_events, random_state=0)
                                                 for level in level_combinations],
                                                axis=0)

    # anova selection_dict create and save
    anova_selection_dict = {'df': df,
                            'selection_df': selection_df,
                            'counts': counts,
                            'num_events': num_events,
                            'valid': valid}
    md.np_saver(anova_selection_dict, target_filename)

main()
