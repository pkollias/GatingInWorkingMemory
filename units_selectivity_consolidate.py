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

        # anova_time_results file load parameters
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

        physiology_df = physiology
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


import matplotlib.pyplot as plt

def cumul_distr_from_series(series):
    cf_count = scipy.stats.cumfreq(np.array(series),
                                   numbins=1000,
                                   defaultreallimits=(0, 1)).cumcount
    return (cf_count - min(cf_count)) / (max(cf_count) - min(cf_count))

df = pd.concat([units, physiology_df], axis=1)

for area in ['PFC', 'Stri', 'IT']:

    df_area = df.loc[df['Area'].eq(area)]

    # Gating modulation scatter plot
    color_list = ['#d817a5', '#48b213', '#727272']
    line_list = ['PFC', 'Stri', 'IT']
    color = dict(zip(line_list, color_list))


    ax1 = df_area.plot.scatter(x='GatingModulation', y='PreDistModulation', c=color[area])
    lims = (-.5, 5.5)
    ax1.set_ylim(lims)
    ax1.set_xlim(lims)
    ax1.plot(list(lims), list(lims), color='k')
    ax1.plot(list(lims), [1, 1], color='k', ls=':')
    ax1.plot([1, 1], list(lims), color='k', ls=':')
    ax1.set_xlabel('Gating Modulation ($FR_{gating}/FR_{fixation}$)')
    ax1.set_ylabel('PreDist Modulation ($FR_{predist}/FR_{fixation}$)')
    ax1.set_title(area)
    plt.savefig('Modulation_{0:s}.png'.format(area), format='png')
    plt.close()


# Depth of selectivity cum hist
for period in ['Overall', 'Distractor']:

    fig = plt.figure()

    for area in ['PFC', 'Stri', 'IT']:

        period_str = 'DepthOfSelectivityIndex_{0:s}'.format(period)
        df_drop = df_area.drop(df_area.loc[df_area[period_str].gt(1)].index)
        # line = df_drop[period_str].dropna().sort_values().cumsum() / df_drop[period_str].sum()
        ax2.plot(np.linspace(0, 1, 1000), cumul_distr_from_series(df_drop[period_str].dropna()), color=color[area])
        # ax2 = df_drop[period_str].hist(cumulative=True, density=1, bins=2000, color=color[area])

    ax2.set_xlabel('Depth of Stimulus Selectivity Index')
    ax2.set_ylabel('Cumulative Distribution of Population')
    ax2.set_title('{0:s} - {1:s}'.format(area, period))
    plt.savefig('DoSSI_{0:s}_{1:s}.png'.format(area, period), format='png')
    plt.close()

