from aov_stats import *
import matplotlib.pyplot as plt
from versioning import *
from matplotlib.patches import Rectangle
from matplotlib.offsetbox import AnchoredText


def selectivity_plot(subject_list, version_fr, x_cluster_label, x_cluster_label_list, area_list, ylim_arg, title, crop=(-50, 1000)):

    # version parameters
    version_aov = 'GeneralizedGating'
    x_a = 'GatingCondSpecialized'
    x_b = 'StageStimSpecialized'
    x_ab = interaction_term(x_a, x_b)
    x_factors = [x_a, x_b, x_ab]

    # subject_list = ['Oscar' 'Gonzo']
    # version_fr = 'ConcatFactor'
    # x_cluster_label = x_a
    # x_cluster_label_list = [(x_a, True), (x_ab, False)]
    # area_list = ['PFC']
    # ylim_arg = None
    # title = ''

    v_fr_params = version_fr_params(version_fr)
    t_start = v_fr_params['t_start']
    t_end = v_fr_params['t_end']
    timebin = v_fr_params['timebin']
    timestep = v_fr_params['timestep']
    timebin_interval = TimebinInterval(timebin, timestep, t_start, t_end)

    selection_plot_list = ['All']
    selection_list = ['All']

    # init
    md = MetaData()
    db = md.db_base_loader(['sessions', 'units'])
    sessions, units = db['sessions'], db['units']
    src_filename = md.proc_dest_path(path.join('BehavioralUnits', 'Anova', version_aov,
                                               behunit_params_str(version_fr, timebin, timestep, t_start, t_end), 'consolidate'),
                                     'physiology_dict.pkl')

    physiology_dict = md.np_loader(src_filename)

    # plotting parameters
    color_list = ['#d817a5', '#48b213', '#727272']
    shade_color_list = ['#f230be', '#5fce27', '#9b9b9b']
    linewidth_list = [2.5, 2.5, 2.5]
    linestyle_list = [(0, ()), (0, ()), (0, ())]
    line_list = ['PFC', 'Stri', 'IT']
    color = dict(zip(line_list, color_list))
    shade_color = dict(zip(line_list, shade_color_list))
    linewidth = dict(zip(line_list, linewidth_list))
    linestyle = dict(zip(line_list, linestyle_list))

    def apply_label_constraint_to_unit(unit, x_cl_list):
        return all([unit[x_label]['valid'] and bool(unit[x_label]['clusters']) == val for x_label, val in x_cl_list])

    ### TODO: polish later. Rushed
    # filter subject units
    units = units.loc[sessions.loc[sessions['Subject'].isin(subject_list)].index]
    # filter only single units of area of interest
    units.drop(units.loc[units['UnitNum'].eq(0) | units['RatingCode'].eq(7) |
                         ~units['Area'].isin(area_list)].index, inplace=True)
    # convert units index into tuple
    units['tuple_index'] = units.apply(lambda row: row.name, axis=1)
    units.set_index('tuple_index', drop=True, inplace=True)
    ### TODO: improve all that
    filtered_ind = set(units.index)
    valid_ind = set.intersection(*[set([unit_ind
                                        for unit_ind, unit in selection_dict.items()
                                        if unit[x_cluster_label]['valid']])
                                   for selection_dict in physiology_dict.values()])
    cluster_valid_ind = set.intersection(*[set([unit_ind
                                        for unit_ind, unit in selection_dict.items()
                                        if apply_label_constraint_to_unit(unit, x_cluster_label_list)])
                                   for selection_dict in physiology_dict.values()])
    filtered_valid_ind = set.intersection(filtered_ind, valid_ind)

    # calculate selective_interaction_ind from selective_population_ind by running this script with interaction x_ab cluster_label
    # filtered_valid_ind = set(filtered_valid_ind).difference(selective_interaction_ind)

    selective_per_selection_ind = dict([(selection, set([unit_ind for unit_ind, unit in selection_dict.items() if bool(unit[x_cluster_label]['clusters'])]).intersection(cluster_valid_ind)) for selection, selection_dict in physiology_dict.items()])
    selective_population_ind = set.union(*[selection_ind for selection_ind in selective_per_selection_ind.values()])
    units_per_area_ind = dict([(area, set(units.loc[units['Area'].eq(area)].index).intersection(filtered_valid_ind)) for area in area_list])



    pev_omega_sq = {}
    # create subdir, area combination collections of omega_sq pev
    for selection in selection_list:
        for area in area_list:

            unit_ind_list = set.intersection(units_per_area_ind[area], selective_population_ind)
            # unit_ind_list = units_per_area_ind[area]

            # get physiology slice
            pev_omega_sq[(selection, area)] = {'array': [physiology_dict[selection][unit_index][x_cluster_label]['zscores']
                                                         for unit_index in unit_ind_list]}
            pev_omega_sq[(selection, area)]['mean'] = np.nanmean(pev_omega_sq[(selection, area)]['array'], axis=0)
            pev_omega_sq[(selection, area)]['std'] = np.nanstd(pev_omega_sq[(selection, area)]['array'], axis=0)
            pev_omega_sq[(selection, area)]['n'] = np.array(pd.DataFrame(pev_omega_sq[(selection, area)]['array']).notna().sum())
            pev_omega_sq[(selection, area)]['se'] = pev_omega_sq[(selection, area)]['std']/np.sqrt(pev_omega_sq[(selection, area)]['n'])

    ############
    ############
    # plot
    smooth_filter = signal.windows.gaussian(2, .5)
    smoother = SignalSmoothing(signal.correlate, smooth_filter).smoothen # lambda x: x # SignalSmoothing(signal.correlate, signal.windows.gaussian(4, 1)).smoothen

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    legend = []

    for selection in selection_plot_list:
        for area in area_list:

            plot_dim = area

            clr = color[plot_dim]
            shd_clr = shade_color[plot_dim]
            ls = linestyle[plot_dim]
            lw = linewidth[plot_dim]

            significant_n = len(set.intersection(selective_per_selection_ind[selection], units_per_area_ind[area]))
            selective_group_n = len(set.intersection(selective_population_ind, units_per_area_ind[area]))
            area_n = len(units_per_area_ind[area])

            tb_columns = timebin_interval.split_to_bins_offset()
            t_slice = slice(tb_columns.index(crop[0]), tb_columns.index(crop[1]) + 1)
            tb_columns = tb_columns[t_slice]
            sig_mean = smoother(pev_omega_sq[(selection, area)]['mean'])[t_slice]
            sig_se = smoother(pev_omega_sq[(selection, area)]['se'])[t_slice]

            ax.plot(tb_columns, sig_mean, color=clr, ls=ls, lw=lw)
            ax.fill_between(tb_columns, sig_mean - sig_se, sig_mean + sig_se, alpha=0.2,
                            edgecolor=shd_clr,
                            facecolor=shd_clr)

            subdir_str = '' if selection == 'All' else selection
            legend.append('{0:s} {1:s} (N={2:d}, {3:d}%)'.format(area, subdir_str, significant_n, round(100*(float(significant_n)/float(area_n)))))

    ax.legend(legend, prop={'size': 7})
    ylims = ax.get_ylim() if not ylim_arg else ylim_arg
    plt.ylim(ylims[0], ylims[1])
    xlims = ax.get_xlim()
    ax.add_patch(Rectangle((0, ylims[0]), 385, ylims[1]-ylims[0], color='k', alpha=0.1))
    ax.hlines(0, xlims[0], xlims[1], '#999999', linestyles='dashed')

    # title_str = '{0:s}: {1:s}'.format('_'.join(subject_list), title)
    title_str = '{1:s}'.format('_'.join(subject_list), title)
    ax.set_title(title_str)
    ax.set_xlabel('event time (ms) relative to sample onset')
    ax.set_ylabel('z-scored $\omega^2$')
    ax.autoscale(enable=True, axis='x', tight=True)

    return ax


def scatter_plot(subject_list, area, version_fr, selective):

    v_fr_params = version_fr_params(version_fr)
    t_start = v_fr_params['t_start']
    t_end = v_fr_params['t_end']
    timebin = v_fr_params['timebin']
    timestep = v_fr_params['timestep']

    # version parameters
    version_aov = 'GeneralizedGating'
    x_a = 'GatingCondSpecialized'
    x_b = 'StageStimSpecialized'
    x_ab = interaction_term(x_a, x_b)
    x_factors = [x_a, x_b, x_ab]
    selection = 'All'
    x_x = x_ab
    y_x = x_a

    # init
    rounding_decimal = 10
    anova_round = lambda x: round(x, rounding_decimal)

    md = MetaData()
    db = md.db_base_loader(['sessions', 'units'])
    sessions, units = db['sessions'], db['units']


    units = units.loc[sessions.loc[sessions['Subject'].isin(subject_list)].index]
    # filter only single units of area of interest
    units.drop(units.loc[units['UnitNum'].eq(0) | units['RatingCode'].eq(7) |
                         ~units['Area'].eq(area)].index, inplace=True)
    # convert units index into tuple
    units['tuple_index'] = units.apply(lambda row: row.name, axis=1)
    units.set_index('tuple_index', drop=True, inplace=True)
    filtered_ind = set(units.index)


    zscore_unit_results = {}
    for unit_ind in filtered_ind:

        sess, channum, unitnum = unit_ind
        src_filename = md.proc_dest_path(path.join('BehavioralUnits', 'Anova', version_aov,
                                                   behunit_params_str(version_fr, timebin, timestep, t_start, t_end),
                                                   'time_anova', selection),
                                         'time_anova_results_{0:s}_chan{1:03d}_unit{2:03d}.pkl'.
                                         format(sess, channum, unitnum))
        unit_time_results = md.np_loader(src_filename)

        if not bool(unit_time_results):
            print(unit_ind, 'Fail')
            continue
        else:
            print(unit_ind, 'Pass')

        zscore_unit_results[unit_ind] = {}
        for x_i in [x_x, y_x]:

            shuffles_series = unit_time_results[y_bin_strnum_to_str(str(t_start))]['results']['shuffles'][x_i]
            omega_sq_distr = shuffles_series.apply(anova_round)
            omega_sq_observed = omega_sq_distr.loc[shuffle_to_name(0)]
            shuffle_values = omega_sq_distr.loc[~omega_sq_distr.index.isin([shuffle_to_name(0)])]
            omega_sq_mean = np.mean(shuffle_values)
            omega_sq_std = np.std(shuffle_values)
            zscore = anova_round((omega_sq_observed - omega_sq_mean) / omega_sq_std) if omega_sq_std > 0 else 0
            zscore_unit_results[unit_ind][x_i] = zscore


    x, y = zip(*[(val[x_x], val[y_x]) for val in zscore_unit_results.values()])
    x_thr = 1.96 #np.percentile(x, 95)
    y_thr = 1.96 #np.percentile(y, 95)

    q = [np.ravel_multi_index([int(x_i > x_thr), int(y_i > y_thr)], (2, 2)) for x_i, y_i in zip(x, y)]
    # 1 | 3
    # ----- quadrants
    # 0 | 2
    color = ['#a85e32', '#3f9c33', '#9c3589', '#3768a3']
    c = [color[q_i] for q_i in q]
    df = pd.DataFrame(zip(x, y, q, c), columns=['x', 'y', 'q', 'c'])




    # fig, axs = plt.subplots(2, 2, figsize=(15, 15))

    # ax = axs[0, 0]
    fig, ax = plt.subplots(1, 1)
    plt.sca(ax)
    plt.scatter(x, y, c=c)
    xlim = plt.xlim()
    ylim = plt.ylim()
    plt.hlines(y_thr, xlim[0], xlim[1], '#999999', linestyles='dashed')
    plt.vlines(x_thr, ylim[0], ylim[1], '#999999', linestyles='dashed')
    plt.plot([max(xlim[0], ylim[0]), min(xlim[1], ylim[1])], [max(xlim[0], ylim[0]), min(xlim[1], ylim[1])], color='#999999')
    loc_list = ['lower left', 'upper left', 'lower right', 'upper right']
    qhist = np.histogram(q, [0, 1, 2, 3, 4], density=True)[0]
    for at_i in range(4):
        at = AnchoredText('{0:.2f}'.format(qhist[at_i]), loc=loc_list[at_i], prop=dict(size=9, color=color[at_i], weight='bold'), frameon=False)
        ax.add_artist(at)
    # plt.hlines(ylim[0], xlim[0], xlim[1], '#a12d3a', linewidth=3)
    # plt.hlines(ylim[1], xlim[0], xlim[1], '#22b5a6', linewidth=3)
    # plt.vlines(xlim[0], ylim[0], ylim[1], '#77942c', linewidth=3)
    # plt.vlines(xlim[1], ylim[0], ylim[1], '#582c91', linewidth=3)
    ax.set_xlabel('zscored $\omega^2_{interaction}$')
    ax.set_ylabel('zscored $\omega^2_{Gating}$')


    # ax = axs[1, 0]
    # plt.sca(ax)
    # histplot(data=df.loc[df['q'].isin([1, 3])], x='x', stat='probability', element='step', binwidth=5, color='#22b5a6', shrink=0.8)
    # histplot(data=df.loc[df['q'].isin([0, 2])], x='x', stat='probability', element='step', binwidth=5, color='#a12d3a', shrink=0.8)
    # ax.set_xlabel('zscored $\omega^2_{interaction}$')
    # ax.legend(['<.05 $\omega^2_{Gating}$', 'No Gating'])
    #
    # ax = axs[0, 1]
    # plt.sca(ax)
    # histplot(data=df.loc[df['q'].isin([2, 3])], y='y', stat='probability', element='step', binwidth=5, color='#582c91', shrink=0.8)
    # histplot(data=df.loc[df['q'].isin([0, 1])], y='y', stat='probability', element='step', binwidth=5, color='#77942c', shrink=0.8)
    # ax.set_ylabel('zscored $\omega^2_{Gating}$')
    # ax.legend(['<.05 $\omega^2_{interaction}$', 'No interaction'])
    #
    # ax = axs[1, 1]
    # plt.sca(ax)
    # ax.set_visible(False)

    fig.suptitle(area, fontsize=12)
    fig.savefig('GeneralizedGating_{0:s}_{1:s}_Hist.png'.format(area, version_fr), format='png')

