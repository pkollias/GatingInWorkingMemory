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

    columns = ['GatingModulation', 'PreDistModulation', 'GatingSelectivityIndex',
               'StimulusSelectivityIndex_Overall', 'DepthOfSelectivityIndex_Overall',
               'StimulusSelectivityIndex_Distractor', 'DepthOfSelectivityIndex_Distractor']


    physiology_dict = {}

    for version_fr in ['WindowSampleShift', 'WindowDelayShift']:

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

            unit_selectivity_entry = [np.nan for _ in range(7)]
            if bool(unit_selectivity_index):
                unit_selectivity_entry = [unit_selectivity_index['modulation_ratio']['Gating'],
                                          unit_selectivity_index['modulation_ratio']['PreDist'],
                                          unit_selectivity_index['modulation_ratio']['GatingSelectivityIndex'],
                                          unit_selectivity_index['selectivity_index']['Overall']['StimulusSelectivityIndex'],
                                          unit_selectivity_index['selectivity_index']['Overall']['DepthOfSelectivityIndex'],
                                          unit_selectivity_index['selectivity_index']['Distractor']['StimulusSelectivityIndex'],
                                          unit_selectivity_index['selectivity_index']['Distractor']['DepthOfSelectivityIndex']]

            physiology_dict[unit_ind] = unit_selectivity_entry

        physiology_df = pd.DataFrame.from_dict(physiology_dict, orient='index', columns=columns)
        physiology_df.index = pd.MultiIndex.from_tuples(physiology_df.index)
        physiology_df.index.names = physiology.index.names
        physiology_df = pd.concat([physiology, physiology_df], axis=1)

        md.np_saver(physiology_df, target_filename)



main()
