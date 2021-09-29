from rec import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
from matplotlib.patches import Circle
from scipy.stats import binom_test
from scipy.stats import ttest_1samp
from rec_utils import *


def selectivity_plot(subject_list, version_aov, version_fr, selection_plot_list, selection_list, area_list,
                     ylim_arg, selective_across, title, legend_show, legend_cond_show, percent_show,
                     cluster_corrected, plot_nonzero, ax, alpha=(1, 0.2), crop=(-50, 1000), ax2_y=None):

    # subject_list = ['Oscar', 'Gonzo']
    # version_aov = 'PresentedStimulus'
    # version_fr = 'ConcatFactor'
    # selection_list = ['PreDist', 'Gating', 'PostDist', 'Target']
    # selection_plot_list = selection_list
    # area_list = ['PFC']
    # ylim_arg = None
    # title = ''

    ax_passed = not ax is None

    # version parameters
    v_aov_params = anova_version_aov_params(version_aov, version_fr)
    x_factors = v_aov_params['x_factors']
    x_cluster_label = x_factors[0]

    v_fr_params = version_fr_params(version_fr)
    t_start = v_fr_params['t_start']
    t_end = v_fr_params['t_end']
    timebin = v_fr_params['timebin']
    timestep = v_fr_params['timestep']

    # init
    md = MetaData()
    db = md.db_base_loader(['sessions', 'units'])
    sessions, units = db['sessions'], db['units']
    file_str = 'consolidate' if cluster_corrected else 'summarize'
    # src_filename = md.proc_dest_path(path.join('BehavioralUnits', 'Anova', version_aov,
    #                                            behunit_params_str(version_fr, timebin, timestep, t_start, t_end), file_str),
    #                                  'physiology_dict.pkl')


    anova = Anova(DataBase([]), {'aov': version_aov, 'fr': version_fr})
    anova.load_physiology_dict()
    physiology_dict = anova.physiology_dict


    if not bool(selection_list) and not bool(selection_plot_list):
        selection_list = v_aov_params['selection_dict']['list']
        selection_plot_list = v_aov_params['selection_dict']['list']

    # plotting parameters
    if selection_list == ['All']:
        color_list = ['#d817a5', '#48b213', '#727272']
        shade_color_list = ['#f230be', '#5fce27', '#9b9b9b']
        linewidth_list = [2.5, 2.5, 2.5]
        linestyle_list = [(0, ()), (0, ()), (0, ())]
        line_list = ['PFC', 'Stri', 'IT']
        color = dict(zip(line_list, color_list))
        shade_color = dict(zip(line_list, shade_color_list))
        linewidth = dict(zip(line_list, linewidth_list))
        linestyle = dict(zip(line_list, linestyle_list))
    else:
        color_list = ['#317275', '#dbb356', '#3c3a7a', '#e08e4f', '#808080', '#ffb778', '#e08e4f', '#e08e4f', '#e08e4f',
                      '#bd9842', '#282761', '#bd6b2d', '#ab3a0e', '#9e9e9e']
        shade_color_list = ['#49898c', '#e8c36d', '#504dab', '#f0a869', '#a8a8a8', '#eda768', '#f0a869', '#f0a869', '#f0a869',
                            '#d1a94b', '#323073', '#db8340', '#c74816', '#a8a8a8']
        linewidth_list = [2, 2, 2.5, 2, 2, 2.5, 2, 2, 2,
                          2, 2, 2, 2, 2]
        # linestyle_list = [(0, ()), (0, (2, 1, 1, 1)), (0, ()), (0, (5, 4)), (0, (5, 1)), (0, ()), (0, (1, 1)), (0, (2, 1)), (0, (3, 1)),
        #                   (0, ()), (0, ()), (0, ()), (0, ()), (0, ())]
        linestyle_list = [(0, ()), (0, ()), (0, ()), (0, ()), (0, ()), (0, ()), (0, ()), (0, ()), (0, ()),
                          (0, ()), (0, ()), (0, ()), (0, ()), (0, ())]
        line_list = ['Cue', 'PreDist', 'Gating', 'PostDist', 'Target', 'PreDist_PostDist', 'PostDist1', 'PostDist2', 'PostDist3',
                     'PreDist_From_To_PreDist', 'PreDist_From_To_Gating', 'Gating_From_To_PostDist', 'PostDist_From_To_PostDist', 'PostDist_From_To_Target']
        color = dict(zip(line_list, color_list))
        shade_color = dict(zip(line_list, shade_color_list))
        linewidth = dict(zip(line_list, linewidth_list))
        linestyle = dict(zip(line_list, linestyle_list))

    timebin_interval = TimebinInterval(timebin, timestep, t_start, t_end)




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
    filtered_inds = set(units.index)
    valid_inds = set.intersection(*[set([unit_ind for unit_ind, unit in selection_dict.items() if unit[x_cluster_label]['valid']]) for selection_dict in physiology_dict.values()])
    filtered_valid_inds = set.intersection(filtered_inds, valid_inds)

    if cluster_corrected:
        selective_per_selection_inds = dict([(selection, set([unit_ind for unit_ind, unit in selection_dict.items() if bool(unit[x_cluster_label]['clusters'])]).intersection(filtered_valid_inds)) for selection, selection_dict in physiology_dict.items() if selection in selection_list])
        selective_population_inds = set.union(*[selection_ind for selection_ind in selective_per_selection_inds.values()])
    units_per_area_inds = dict([(area, set(units.loc[units['Area'].eq(area)].index).intersection(filtered_valid_inds)) for area in area_list])

    # get p values of binomial test for selective units for units x timebins selectivity array
    def get_bins_pvals_binom(selectivity_array):
        return [binom_test(list(pd.Series(selectivity_array[:, ii]).dropna().value_counts().reindex([1.0, 0.0]).replace({np.nan: 0})),
                           p=0.05, alternative='greater')
                for ii in range(selectivity_array.shape[1])]

    def get_bins_pvals_ttest_zero(z_score_array):
        return ttest_1samp(np.array(z_score_array), 0, axis=0)[1]


    def pval_to_significance_level(pval):
        if 0.05 < pval < np.inf:
            return 0
        elif 0.01 < pval <= 0.05:
            return 1
        elif 0.001 < pval <= 0.01:
            return 2
        else:
            return 3

    # omega squared
    pev_omega_sq = {}
    # create subdir, area combination collections of omega_sq pev and bins pvals
    for selection in selection_list:
        for area in area_list:

            if cluster_corrected and selective_across:
                unit_ind_list = set.intersection(units_per_area_inds[area], selective_population_inds)
            else:
                unit_ind_list = units_per_area_inds[area]

            # get physiology slice
            pev_omega_sq[(selection, area)] = {'array': [physiology_dict[selection][unit_index][x_cluster_label]['zscores']
                                                         for unit_index in unit_ind_list]}
            pev_omega_sq[(selection, area)]['mean'] = np.nanmean(pev_omega_sq[(selection, area)]['array'], axis=0)
            pev_omega_sq[(selection, area)]['std'] = np.nanstd(pev_omega_sq[(selection, area)]['array'], axis=0)
            pev_omega_sq[(selection, area)]['n'] = np.array(pd.DataFrame(pev_omega_sq[(selection, area)]['array']).notna().sum())
            pev_omega_sq[(selection, area)]['se'] = pev_omega_sq[(selection, area)]['std']/np.sqrt(pev_omega_sq[(selection, area)]['n'])

    # significance
    bins_pvals = {}
    for selection in selection_list:
        for area in area_list:
            if cluster_corrected:
                unit_ind_list = selective_per_selection_inds[selection]

            # get physiology slice
            array = [physiology_dict[selection][unit_index][x_cluster_label]['zscores'] for unit_index in unit_ind_list]
            significance = get_bins_pvals_ttest_zero(array)
            # significance = get_bins_pvals_binom(anova.get_unit_time_selectivity(selection, x_cluster_label, timebin_interval, filtered_inds))
            bins_pvals[(selection, area)] = list(map(pval_to_significance_level, significance))

    ############
    ############
    # plot
    # fig = plt.figure(figsize=(7.8, 4.8))
    if not bool(ax):
        if plot_nonzero:
            fig = plt.figure(figsize=(9.55, 5.35)) if version_fr in ['ConcatFactorExtended'] else plt.figure(figsize=(6.55, 5.35))
            gs = fig.add_gridspec(nrows=40, ncols=1, left=0.15, right=0.7, bottom=0.15, top=0.75, wspace=0.05)
            ax = fig.add_subplot(gs[5:, 0])
            ax2 = fig.add_subplot(gs[:5, 0], sharex=ax)
        else:
            fig = plt.figure(figsize=(8.55, 4.35)) if version_fr in ['ConcatFactorExtended'] else plt.figure(figsize=(5.55, 4.35))
            ax = fig.add_subplot(1, 1, 1)

    legend = []
    if version_fr in ['ConcatFactorExtended']:
        smooth_filter = signal.windows.gaussian(4, 2.3)
    else:
        smooth_filter = signal.windows.gaussian(2, 1)
    smoother = SignalSmoothing(signal.correlate, smooth_filter).smoothen # lambda x: x # SignalSmoothing(signal.correlate, signal.windows.gaussian(4, 1)).smoothen


    # print(ax_passed, type(ax), len(ax))
    if ax_passed and type(ax) == tuple and len(ax) == 2:
        ax2 = ax[1]
        ax = ax[0]

    for selection in selection_plot_list:
        for area in area_list:

            if selection == 'All':
                plot_dim = area
            else:
                plot_dim = selection

            full_tb_columns = timebin_interval.split_to_bins_offset()
            t_slice = slice(full_tb_columns.index(crop[0]), full_tb_columns.index(crop[1]) + 1)
            tb_columns = full_tb_columns[t_slice]
            sig_mean = smoother(pev_omega_sq[(selection, area)]['mean'])[t_slice]
            sig_se = smoother(pev_omega_sq[(selection, area)]['se'])[t_slice]
            clr = color[plot_dim]
            shd_clr = shade_color[plot_dim]
            ls = linestyle[plot_dim]
            lw = linewidth[plot_dim]

            if cluster_corrected:
                significant_n = len(set.intersection(selective_per_selection_inds[selection], units_per_area_inds[area]))
                selective_group_n = len(set.intersection(selective_population_inds, units_per_area_inds[area]))
            area_n = len(units_per_area_inds[area])

            ax.plot(tb_columns, sig_mean, color=clr, ls=ls, lw=lw, alpha=alpha[0])
            ax.fill_between(tb_columns, sig_mean - sig_se, sig_mean + sig_se, alpha=alpha[1],
                            edgecolor=shd_clr,
                            facecolor=shd_clr)

            subdir_str = selection if legend_cond_show else ''
            legend_str = '{0:s}'.format(subdir_str)
            if cluster_corrected and percent_show:
                legend_str = '{1:s} (N={2:d}, {3:d}%)'.format(area, subdir_str, significant_n,
                                                                    round(100 * (float(significant_n) / float(area_n))))
            legend_element = Circle((0, 0), radius=1, edgecolor='w', facecolor=clr, label=legend_str)
            legend.append(legend_element)

    ylims = ax.get_ylim() if not ylim_arg else ylim_arg
    plt.ylim(ylims[0], ylims[1])
    xlims = ax.get_xlim()
    if not ax_passed:
        ax.add_patch(Rectangle((0, ylims[0]), 385, ylims[1]-ylims[0], color='k', alpha=0.1))
        if crop[1] > 950 + 385:
            ax.add_patch(Rectangle((950, ylims[0]), 385, ylims[1] - ylims[0], color='k', alpha=0.1))
        ax.hlines(0, tb_columns[0], tb_columns[-1], '#999999')
    if bool(title):
        title_str = title
    else:
        if cluster_corrected:
            title_str = '{1:s} (N={2:d}/{3:d})'.format('_'.join(subject_list), '_'.join(area_list), selective_group_n, area_n)
        else:
            title_str = '{1:s} (N={2:d})'.format('_'.join(subject_list), '_'.join(area_list), area_n)
    if crop[1] > 950:
        ax.set_xticks(list(np.linspace(0, 600, 3)) + list(950 + np.linspace(0, 600, 3)))
    else:
        ax.set_xticks(list(np.linspace(0, 900, 4)))
    ax.set_xlabel('time since sample onset ($ms$)')
    ax.set_ylabel('z-scored $\omega^2$')

    ###

    if plot_nonzero:
        plt.sca(ax2)
        ax2.set_ylim((0 - .5, (len(selection_list) if not ax2_y else (ax2_y[1] + 1)) - 1 + .5))
        plt.axis('off')
        for sel_ii, selection in enumerate(selection_plot_list):
            for area in area_list:
                x = np.array(tb_columns)
                y = np.ones_like(x) * sel_ii if not ax2_y else np.ones_like(x) * ax2_y[0]
                lwidths = np.array(bins_pvals[(selection, area)])[t_slice] * 1
                points = np.array([x, y]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                segments = segments # + np.tile([[1, 0], [-1, 0]], (40, 1, 1)) * 0.3
                lc = LineCollection(segments, linewidths=lwidths, color=color[selection], alpha=alpha[0])
                ax2.add_collection(lc)

    ###

    if legend_show:
        ax.legend(handles=legend, loc='upper right', bbox_to_anchor=(1, 1), frameon=False, prop={'size': 9})

    # ax.set_position([0.18, 0.1, 0.54, 0.8])
    if legend_show:
        ax.get_legend().set_bbox_to_anchor((1, .5))
        ax.get_legend()._set_loc(6)
    # ax.autoscale(enable=True, axis='x', tight=True)
    # plt.tight_layout()
    plt.sca(ax)
    plt.box(False)

    if plot_nonzero:
        ax2.set_title(title_str)
        return (ax, ax2)
    else:
        ax.set_title(title_str)
        return ax
