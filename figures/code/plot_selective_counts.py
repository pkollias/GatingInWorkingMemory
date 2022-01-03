from rec import *
import matplotlib.pyplot as plt
from rec_utils import *
import math
from matplotlib.lines import Line2D

def selectivity_count_plot(subject_list, version_aov, version_fr, selection_list, area_list):

    # subject_list = ['Oscar', 'Gonzo']
    # version_aov = 'PresentedStimulus'
    # version_fr = 'ConcatFactor'
    # selection_list = ['PreDist', 'Gating', 'PostDist', 'Target']
    # area_list = ['PFC']
    # ylim_arg = None
    # title = ''

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

    anova = Anova(DataBase([]), {'aov': version_aov, 'fr': version_fr})
    anova.load_physiology_dict(cropped=True)
    physiology_dict = anova.physiology_dict

    if not bool(selection_list):
        selection_list = v_aov_params['selection_dict']['list']

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

    selective_per_selection_inds = dict([(selection, set([unit_ind for unit_ind, unit in selection_dict.items() if bool(unit[x_cluster_label]['clusters'])]).intersection(filtered_valid_inds)) for selection, selection_dict in physiology_dict.items() if selection in selection_list])
    selective_population_inds = set.union(*[selection_ind for selection_ind in selective_per_selection_inds.values()])
    units_per_area_inds = dict([(area, set(units.loc[units['Area'].eq(area)].index).intersection(filtered_valid_inds)) for area in area_list])

    fig = plt.figure(figsize=(1.47, 2.91))
    ax = fig.add_subplot(1, 1, 1)

    selective_n = {}
    area_n = {}
    for area in area_list:

        selective_n = {}
        area_n = len(units_per_area_inds[area])

        for selection in selection_list:
            selective_n[selection] = len(set.intersection(selective_per_selection_inds[selection], units_per_area_inds[area])) / area_n * 100

        selective_n['All'] = len(set.intersection(*[set.intersection(selective_per_selection_inds[selection], units_per_area_inds[area])
                                                          for selection in selection_list])) / area_n * 100
        selective_n['Any'] = len(set.union(*[set.intersection(selective_per_selection_inds[selection], units_per_area_inds[area])
                                                   for selection in selection_list])) / area_n * 100

    print([el  * 100 * area_n for el in selective_n.values()])
    plt.bar(range(len(selective_n)), selective_n.values())
    ax.set_title(area)
    ax.set_ylim(0, 35)
    ax.set_ylabel('Selective Units (%)')
    plt.box(False)

    return fig, ax


def selectivity_counts(subject_list, version_aov, version_fr, area, selection_list, selection_mode='Single'):

    # subject_list = ['Oscar', 'Gonzo']
    # version_aov = 'PresentedStimulus'
    # version_fr = 'ConcatFactor'
    # selection_list = ['PreDist', 'Gating', 'PostDist', 'Target']
    # area_list = ['PFC']
    # ylim_arg = None
    # title = ''

    area_list = [area]

    # version parameters
    v_aov_params = anova_version_aov_params(version_aov, version_fr)
    x_factors = v_aov_params['x_factors']
    x_cluster_label = x_factors[0]

    # init
    md = MetaData()
    db = md.db_base_loader(['sessions', 'units'])
    sessions, units = db['sessions'], db['units']

    anova = Anova(DataBase([]), {'aov': version_aov, 'fr': version_fr})
    anova.load_physiology_dict(cropped=True)
    physiology_dict = anova.physiology_dict

    if not bool(selection_list):
        selection_list = v_aov_params['selection_dict']['list']


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

    selective_per_selection_inds = dict([(selection, set([unit_ind for unit_ind, unit in selection_dict.items() if bool(unit[x_cluster_label]['clusters'])]).intersection(filtered_valid_inds)) for selection, selection_dict in physiology_dict.items() if selection in selection_list])
    selective_population_inds = set.union(*[selection_ind for selection_ind in selective_per_selection_inds.values()])
    units_per_area_inds = dict([(area, set(units.loc[units['Area'].eq(area)].index).intersection(filtered_valid_inds)) for area in area_list])


    selective_n = {}
    area_n = len(units_per_area_inds[area])

    for selection in selection_list:
        selective_n[selection] = len(set.intersection(selective_per_selection_inds[selection], units_per_area_inds[area])) / area_n * 100

    selective_n['All'] = len(set.intersection(*[set.intersection(selective_per_selection_inds[selection], units_per_area_inds[area])
                                                      for selection in selection_list])) / area_n * 100
    selective_n['Any'] = len(set.union(*[set.intersection(selective_per_selection_inds[selection], units_per_area_inds[area])
                                               for selection in selection_list])) / area_n * 100

    key = selection_list[0] if selection_mode == 'Single' else selection_mode

    return round(selective_n[key] * area_n / 100), area_n


# from statsmodels.stats.multitest import multipletests
# from scipy.stats import binom_test
# counts1 = [selectivity_counts(['Gonzo', 'Oscar'], 'PresentedStimulus', 'ConcatFactor', area, [stage])
#           for area, stage in product(['PFC', 'Stri', 'IT'], ['PreDist', 'Gating', 'PostDist', 'Target'])]
# counts2 = [selectivity_counts(['Gonzo', 'Oscar'], 'GatedStimulus', 'ConcatFactor', area, [stage])
#           for area, stage in product(['PFC', 'Stri', 'IT'], ['Cue', 'PreDist', 'PostDist'])]
# tests = [binom_test(*count, .05, alternative='greater') for count in counts1 + counts2]
# multipletests(tests, method='fdr_bh')



def selectivity_count_intersectional_plot(subject_list, version_aov, version_fr, selection_list, area_list):

    # subject_list = ['Oscar', 'Gonzo']
    # version_aov = 'PresentedStimulus'
    # version_fr = 'ConcatFactor'
    # selection_list = ['PreDist', 'Gating', 'PostDist', 'Target']
    # area_list = ['PFC']
    # ylim_arg = None
    # title = ''

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

    anova = Anova(DataBase([]), {'aov': version_aov, 'fr': version_fr})
    anova.load_physiology_dict(cropped=True)
    physiology_dict = anova.physiology_dict

    if not bool(selection_list):
        selection_list = v_aov_params['selection_dict']['list']

    # plotting parameters

    units = units.loc[sessions.loc[sessions['Subject'].isin(subject_list)].index]
    # filter only single units of area of interest
    units.drop(units.loc[units['UnitNum'].eq(0) | units['RatingCode'].eq(7) |
                         ~units['Area'].isin(area_list)].index, inplace=True)
    # convert units index into tuple
    units['tuple_index'] = units.apply(lambda row: row.name, axis=1)
    units.set_index('tuple_index', drop=True, inplace=True)
    filtered_inds = set(units.index)
    valid_inds = set.intersection(*[set([unit_ind for unit_ind, unit in selection_dict.items() if unit[x_cluster_label]['valid']]) for selection_dict in physiology_dict.values()])
    filtered_valid_inds = set.intersection(filtered_inds, valid_inds)

    selective_per_selection_inds = dict([(selection, set([unit_ind for unit_ind, unit in selection_dict.items() if bool(unit[x_cluster_label]['clusters'])]).intersection(filtered_valid_inds)) for selection, selection_dict in physiology_dict.items() if selection in selection_list])
    selective_population_inds = set.union(*[selection_ind for selection_ind in selective_per_selection_inds.values()])
    units_per_area_inds = dict([(area, set(units.loc[units['Area'].eq(area)].index).intersection(filtered_valid_inds)) for area in area_list])

    ax = {}
    for area in area_list:

        selective_n = {}
        selective_inds = {}
        area_n = len(units_per_area_inds[area])

        for selection in selection_list:
            selective_inds[selection] = set.intersection(selective_per_selection_inds[selection], units_per_area_inds[area])
            selective_n[selection] = len(selective_inds[selection])

        ax[area] = intersectional_pie([units_per_area_inds[area], selective_inds['PreDist'], selective_inds['Gating'], selective_inds['PostDist']],
                                      legend_show=True if area == 'PFC' else False, title_str=area)

    return ax


def intersectional_pie(set_list, title_str='', legend_show=False, R=1):

    pi = math.pi
    # 0:abc, 1:ab, 2:ac, 3:bc, 4:a, 5:b, 6:c, 7:null
    color_list = ['#eeeeee', '#5a914a', '#9c7038', '#754085', '#dbb356', '#3c3a7a', '#e08e4f', '#404040']
    label_list = ['All', 'PreDist-Gating', 'PreDist-PostDist', 'Gating-PostDist', 'PreDist', 'Gating', 'PostDist', 'None']
    color_xy = itemgetter(4, 2, 6, 3, 5, 1)(color_list)

    omega, alpha, beta, gamma = set_list
    s_abc = set.intersection(alpha, beta).intersection(gamma)
    s_ab = set.intersection(alpha, beta).difference(s_abc)
    s_ac = set.intersection(alpha, gamma).difference(s_abc)
    s_bc = set.intersection(beta, gamma).difference(s_abc)
    s_a = alpha.difference(beta).difference(gamma)
    s_b = beta.difference(alpha).difference(gamma)
    s_c = gamma.difference(alpha).difference(beta)
    s_null = omega.difference(alpha).difference(beta).difference(gamma)

    a_omega = len(omega)
    a_abc = len(s_abc) * 5 / a_omega
    a_ab = len(s_ab) * 5 / a_omega
    a_ac = len(s_ac) * 5 / a_omega
    a_bc = len(s_bc) * 5 / a_omega
    a_a = len(s_a) * 5 / a_omega
    a_b = len(s_b) * 5 / a_omega
    a_c = len(s_c) * 5 / a_omega
    a_null = len(s_null) * 5 / a_omega

    a_xy = a_a + a_b + a_c + a_ab + a_ac + a_bc
    a_omega = len(omega)

    r_abc = (a_abc / pi) ** (1/2)
    r_xy = ((a_abc + a_xy) / pi) ** (1/2) - r_abc
    r_null = ((a_abc + a_xy + a_null) / pi) ** (1/2) - r_abc - r_xy
    R = r_abc + r_xy + r_null

    fig, ax = plt.subplots()
    ax.pie([1], radius=r_abc, colors=[color_list[0]], wedgeprops=dict(width=r_abc, edgecolor=(1, 1, 1, 0), linewidth=.33))
    pie_wedge_collection = ax.pie([a_a, a_ac, a_c, a_bc, a_b, a_ab], normalize=True, radius=r_xy + r_abc, colors=list(color_xy),
                                  wedgeprops=dict(width=r_xy, edgecolor=(1, 1, 1, 0), linewidth=.33))
    ax.pie([1], radius=R, colors=[color_list[-1]], wedgeprops=dict(width=r_null, edgecolor=(1, 1, 1, 0), linewidth=.33))


    ax.set(aspect='equal')
    legend_elements = [Line2D([0], [0], marker='o', color=(1, 1, 1, 0), label=l, markerfacecolor=c, markeredgecolor=(1, 1, 1, 0), markersize=12)
                       for c, l in zip(color_list, label_list)]
    if legend_show:
        ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.35, 0.9), frameon=False, prop={'size': 9})
    plt.show()
    ax.set_title(title_str)

    return ax

