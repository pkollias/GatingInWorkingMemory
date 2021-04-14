import sys
from rec import *
from aov_stats import *


def main():

    # parameters
    args = sys.argv
    version_aov = 'GeneralizedGating'


    # init tables and slices
    md = MetaData()
    db = md.db_base_loader(['units', 'physiology'])
    units, physiology = db['units'], db['physiology']

    columns = ['GatingModulation_All', 'PreDistModulation_All', 'GatingSelectivityIndex_All',
               'GatingModulation_S11', 'PreDistModulation_S11', 'GatingSelectivityIndex_S11',
               'GatingModulation_S12', 'PreDistModulation_S12', 'GatingSelectivityIndex_S12',
               'GatingModulation_S21', 'PreDistModulation_S21', 'GatingSelectivityIndex_S21',
               'GatingModulation_S22', 'PreDistModulation_S22', 'GatingSelectivityIndex_S22',
               'GatingDifference_All', 'PreDistDifference_All', 'GatingContrast_All', 'PreDistContrast_All',
               'GatingDifference_S11', 'PreDistDifference_S11', 'GatingContrast_S11', 'PreDistContrast_S11',
               'GatingDifference_S12', 'PreDistDifference_S12', 'GatingContrast_S12', 'PreDistContrast_S12',
               'GatingDifference_S21', 'PreDistDifference_S21', 'GatingContrast_S21', 'PreDistContrast_S21',
               'GatingDifference_S22', 'PreDistDifference_S22', 'GatingContrast_S22', 'PreDistContrast_S22',
               'StimulusSelectivityIndex_Overall', 'DepthOfSelectivityIndex_Overall',
               'StimulusSelectivityIndex_Distractor', 'DepthOfSelectivityIndex_Distractor',
               'StimulusSelectivityIndex_Gating', 'DepthOfSelectivityIndex_Gating']
    stim_str_levels = ['S11_S12_S21_S22', 'S11', 'S12', 'S21', 'S22']
    gating_mod_str_levels = ['Gating', 'PreDist', 'GatingSelectivityIndex']
    gating_diff_str_levels = ['Gating', 'PreDist', 'GatingContrast', 'PreDistContrast']


    physiology_dict = {}

    for version_fr in ['WindowSampleShift', 'WindowDelayShift']:

        if version_fr == 'WindowDelayShift':
            stim_str_levels = ['S11D_S12D_S21D_S22D', 'S11D', 'S12D', 'S21D', 'S22D']

        v_fr_params = anova_version_fr_params(version_fr)
        t_start = v_fr_params['t_start']
        t_end = v_fr_params['t_end']
        timebin = v_fr_params['timebin']
        timestep = v_fr_params['timestep']
        target_filename = md.proc_dest_path(path.join('BehavioralUnits', 'Anova', version_aov,
                                                      behunit_params_str(version_fr, timebin, timestep, t_start, t_end),
                                                      'selectivity_index'),
                                            'physiology_df.pkl')
        print(target_filename)
        if path.exists(target_filename):
            continue

        physiology_dict = {}

        for u_iloc in range(len(units)):

            unit_entry = units.iloc[u_iloc]
            unit_ind = tuple(unit_entry[md.preproc_imports['units']['index']])
            sess, channum, unitnum = unit_ind

            print(u_iloc, unit_ind)

            src_filename = md.proc_dest_path(path.join('BehavioralUnits', 'Anova', version_aov,
                                                       behunit_params_str(version_fr, timebin, timestep, t_start, t_end),
                                                       'selectivity_index', 'unit_selectivity_index'),
                                             'selectivity_index_{0:s}_chan{1:03d}_unit{2:03d}.pkl'.
                                             format(sess, channum, unitnum))
            # load unit_selectivity_index
            unit_selectivity_index = md.np_loader(src_filename)

            unit_selectivity_entry = [np.nan for _ in range(len(columns))]
            if bool(unit_selectivity_index):
                unit_selectivity_entry = [unit_selectivity_index['modulation_ratio'][s][g]
                                          for s, g
                                          in list(product(stim_str_levels, gating_mod_str_levels))] + \
                                         [unit_selectivity_index['modulation_diff'][s][g]
                                          for s, g
                                          in list(product(stim_str_levels, gating_diff_str_levels))] + \
                                         [unit_selectivity_index['selectivity_index']['Overall']['StimulusSelectivityIndex'],
                                          unit_selectivity_index['selectivity_index']['Overall']['DepthOfSelectivityIndex'],
                                          unit_selectivity_index['selectivity_index']['Distractor']['StimulusSelectivityIndex'],
                                          unit_selectivity_index['selectivity_index']['Distractor']['DepthOfSelectivityIndex'],
                                          unit_selectivity_index['selectivity_index']['Gating']['StimulusSelectivityIndex'],
                                          unit_selectivity_index['selectivity_index']['Gating']['DepthOfSelectivityIndex']]

            physiology_dict[unit_ind] = unit_selectivity_entry

        physiology_df = pd.DataFrame.from_dict(physiology_dict, orient='index', columns=columns)
        physiology_df.index = pd.MultiIndex.from_tuples(physiology_df.index)
        physiology_df.index.names = physiology.index.names
        physiology_df = pd.concat([physiology, physiology_df], axis=1)

        md.np_saver(physiology_df, target_filename)



main()
