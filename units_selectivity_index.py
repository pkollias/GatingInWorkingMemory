import sys
from rec import *
from aov_stats import *


def main():

    # parameters
    args = sys.argv
    u_iloc = int(args[1])
    version_aov = 'GeneralizedGating'
    selection = 'All'

    gating_means_label = 'GatingCondSpecialized'
    gating_label = 'Gating'
    dist_label = 'PreDist'
    stim_means_label = 'StageStimSpecialized'
    inter_means_label = interaction_term(gating_means_label, stim_means_label)




    # init tables and slices
    md = MetaData()
    db = md.db_base_loader(['units'])
    units = db['units']
    unit_entry = units.iloc[u_iloc]
    unit_ind = tuple(unit_entry[md.preproc_imports['units']['index']])
    sess, channum, unitnum = unit_ind


    # fix beh_unit FR file load parameters
    version_fr = 'WindowFix'
    v_fr_params = anova_version_fr_params(version_fr)
    t_start = v_fr_params['t_start']
    t_end = v_fr_params['t_end']
    timebin = v_fr_params['timebin']
    timestep = v_fr_params['timestep']
    src_filename = md.proc_dest_path(path.join('BehavioralUnits', 'FiringRates',
                                               behunit_params_str(version_fr, timebin, timestep, t_start, t_end)),
                                     'behunit_FR_{0:s}_chan{1:03d}_unit{2:03d}.pkl'.format(sess, channum, unitnum))
    # load fix period baseline
    timebin_fr_dict = md.np_loader(src_filename)
    FR_base = np.mean([fr[0] for fr in timebin_fr_dict.values()])


    for version_fr in ['WindowSampleShift', 'WindowDelayShift']:

        if version_fr == 'WindowDelayShift':
            gating_label = gating_label + 'Delay'
            dist_label = dist_label + 'Delay'

        unit_selectivity_index = {}

        # anova_time_results file load parameters
        v_fr_params = anova_version_fr_params(version_fr)
        t_start = v_fr_params['t_start']
        t_end = v_fr_params['t_end']
        timebin = v_fr_params['timebin']
        timestep = v_fr_params['timestep']
        src_filename = md.proc_dest_path(path.join('BehavioralUnits', 'Anova', version_aov,
                                                   behunit_params_str(version_fr, timebin, timestep, t_start, t_end),
                                                   'time_anova', selection),
                                         'time_anova_results_{0:s}_chan{1:03d}_unit{2:03d}.pkl'.
                                         format(sess, channum, unitnum))
        target_filename = md.proc_dest_path(path.join('BehavioralUnits', 'Anova', version_aov,
                                                      behunit_params_str(version_fr, timebin, timestep, t_start, t_end),
                                                      'selectivity_index', 'unit_selectivity_index'),
                                            'selectivity_index_{0:s}_chan{1:03d}_unit{2:03d}.pkl'.
                                            format(sess, channum, unitnum))
        print(target_filename)
        if path.exists(target_filename):
            continue

        # load gating anova results
        unit_time_results = md.np_loader(src_filename)
        if bool(unit_time_results):

            # get means from anova
            means = unit_time_results[y_bin_strnum_to_str(t_start)]['results']['means']
            y_bin_str_start = y_bin_strnum_to_str(t_start)


            # Gating selectivity indices
            FR_gating = means[gating_means_label].loc[gating_label][y_bin_str_start]
            FR_dist = means[gating_means_label].loc[dist_label][y_bin_str_start]
            mod_ratio = lambda fr: fr / FR_base
            unit_selectivity_index['modulation_ratio'] = {'Gating': mod_ratio(FR_gating),
                                                          'PreDist': mod_ratio(FR_dist),
                                                          'GatingSelectivityIndex': mod_ratio(FR_gating - FR_dist)}

            # Stimulus selectivity indices
            unit_selectivity_index['selectivity_index'] = {}
            for stim_fr_period in ['Overall', 'Distractor']:

                # get stimulus means for gating-distractor combined
                if stim_fr_period == 'Overall':
                    FR_stim_i_series = means[inter_means_label].groupby(stim_means_label).mean()[y_bin_str_start]
                # get stimulus means for distractor period
                elif stim_fr_period == 'Distractor':
                    FR_stim_i_series = means[inter_means_label].loc[dist_label][y_bin_str_start]

                StimulusSelectivityIndex = (FR_stim_i_series.max() - FR_stim_i_series.min()) / FR_base
                DepthOfSelectivityIndex = (4 - (FR_stim_i_series / FR_stim_i_series.max()).sum()) / 3

                unit_selectivity_index['selectivity_index'][stim_fr_period] = {'StimulusSelectivityIndex': StimulusSelectivityIndex,
                                                                               'DepthOfSelectivityIndex': DepthOfSelectivityIndex,
                                                                               'FR_stim_i_series': FR_stim_i_series}


        md.np_saver(unit_selectivity_index, target_filename)



main()
