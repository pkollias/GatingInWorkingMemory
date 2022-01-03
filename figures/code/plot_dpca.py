from rec_utils import *
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from itertools import product
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Wedge
from operator import itemgetter
import matplotlib.animation as animation

# args_version = ['factor=StimulusGating', 'fr=ConcatFactor2', 'counts_thr=12',
#                 'area_list=PFC', 'subject=Gonzo_Oscar', 'area=PFC', 'mode=Full', 'mode_seed=0']
# version = job_scheduler(args_version)


def dpca_state_space_plot(version, margin_dim_list, condition_list_slice=False, legend_show=False, azim=-60, elev=30,
                          legend_letters=('1', '2'), parallel=True):

    version_factor = version['factor']
    version_fr = version['fr']
    version_areas = version['area']
    version_subjects = version['subject']
    counts_thr = int(version['counts_thr'])
    # margin_dim_list = [('s', 0), ('s', 1), ('s', 2)]

    # version_parameters
    ver_str = 'dPCA_{0:s}_{1:s}_{2:s}_{3:s}_{4:03d}_{5:s}'.format(version_factor, version_fr,
                                                                        version_areas, version_subjects, counts_thr,
                                                                        margin_dim_list_to_str(margin_dim_list))

    dpca = DemixedPrincipalComponent(DataBase([]), version)
    db, md = dpca.db, dpca.db.md
    target_filename = [dpca.get_path_base(fn, dpca.get_exec_stem()) for fn in ['fbt', 'dpca_obj', 'X_tuple']]
    fbt, dpca_obj, (_, X_trial, _, _) = [md.np_loader(fn) for fn in target_filename]

    # plotting parameters
    (linestyle_column, linestyle), (color_column, color), (alpha_column, alpha), comp = plotting_parameters(version_factor)
    condition_labels = fbt['t'].condition_labels
    condition_list = fbt['t'].df['Condition'].unique()

    tb_columns = fbt['t'].timebin_interval.split_to_bins_offset()
    nearest_t_index = lambda t: min(enumerate(np.abs(np.array(tb_columns)-t)), key=itemgetter(1))[0]
    stim_on_index, delay_on_index, delay_off_index = map(nearest_t_index, [0, 385, 950])
    scatter_color = ['g', 'r', 'k']
    n_dims = len(margin_dim_list)

    # utility functions
    def get_condition_margin_dim_timeseries(fbt, cond, margin, dim):
        if dim < 0:
            series = fbt['t'].timebin_interval.split_to_bins_offset()
        else:
            df = fbt[margin].df
            bool_selection = df['Factor'].eq(dim) & df['Condition'].apply(lambda x: x == cond)
            series = df.loc[bool_selection]['Timeseries'].item()
        return series


    # for every condition get line plotting parameters
    plot_params = {}
    for cond in condition_list: # condition_levels ### TODO

        line_points = tuple([get_condition_margin_dim_timeseries(fbt, cond, margin, dim) for margin, dim in margin_dim_list])
        color_cond = color[cond[condition_labels.index(color_column)]]
        linestyle_cond = linestyle[cond[condition_labels.index(linestyle_column)]]
        alpha_cond = alpha[cond[condition_labels.index(alpha_column)]]
        linewidth_cond = 3  if version_factor == 'PostStimulusRuleStim' and cond[0] == cond[1] \
                            or version_factor == 'GatPostStimulusRuleStim' and cond[0] == cond[1] else \
                         4 if version_factor == 'StimulusGating' else 1.5
        plot_params[cond] = {'line_points': line_points, 'color': color_cond,
                             'linestyle': linestyle_cond, 'linewidth': linewidth_cond, 'alpha': alpha_cond}

    # legend
    legend = []
    # # linestyle
    # linestyle_lines = list(set([(k[condition_labels.index(linestyle_column)], v['linestyle']) for k, v in plot_params.items()]))
    # linestyle_lines = sorted(linestyle_lines, key=lambda x: x[0])
    # for ll in linestyle_lines:
    #     legend.append(Line2D([0], [0], color='#888888', lw=4, ls=ll[1], label='-'.join([legend_letters[0], ll[0]])))
    # # color
    # color_lines = list(set([(k[condition_labels.index(color_column)], v['color']) for k, v in plot_params.items()]))
    # color_lines = sorted(color_lines, key=lambda x: x[0])
    # for cl in color_lines:
    #     legend.append(Line2D([0], [0], color=cl[1], lw=4, label='-'.join([legend_letters[1], cl[0]])))

    cond_list = ['PreDist', 'Gating']
    cond_getter = lambda cond: cond[1]

    # plot
    fig = plt.figure(figsize=(6.5, 5))
    ax = plt.subplot('111') if n_dims < 3 else plt.subplot(111, projection='3d')
    plot_method = lambda ax, n_dims: ax.plot if n_dims < 3 else ax.plot3D
    scatter_method = lambda ax, n_dims: ax.scatter if n_dims < 3 else ax.scatter3D
    task_scatter_points = lambda line_points: [itemgetter(stim_on_index, delay_on_index, delay_off_index)(dim_points) for dim_points in line_points]
    # plot
    for cond, line_params in plot_params.items():
        # if cond_getter(cond) in cond_list:
        plot_method(ax, n_dims)(*line_params['line_points'], color=line_params['color'],
                                linestyle=line_params['linestyle'], linewidth=line_params['linewidth'], alpha=line_params['alpha'])
        scatter_method(ax, n_dims)(*task_scatter_points(line_params['line_points']), color=scatter_color, s=5)


    if n_dims == 3:
        ax.view_init(elev=elev, azim=azim)

    if legend_show:
        leghand = ax.legend(handles=legend, loc='upper right', bbox_to_anchor=(1, 1), frameon=False, prop={'size': 9})
        for legobj in leghand.legendHandles:
            legobj.set_linewidth(3.0)
    label_list = ['Time ($ms$)' if md[1] < 0 else '{0:s}-PC{1:d}'.format(md[0], md[1] + 1) for md in margin_dim_list]
    ax.set_xlabel(label_list[0])
    ax.set_ylabel(label_list[1])
    if n_dims == 3:
        ax.set_zlabel(label_list[2])
    # plt.title(' - '.join([comp[md[0]] + ' PC' + str(md[1]) for md in margin_dim_list]))

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    if n_dims == 3:
        zlim = ax.get_zlim()

    if n_dims == 2:
        if margin_dim_list[0][1] < 0:
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

    ax1 = ax
    plt.savefig(path.join('figures', 'dPCA', '{0:s}.png'.format(ver_str)), format='png')
    # plt.close()



    ####
    # ANIMATION
    ####

    # convert plot_params if all margin_dim_list == 'g'
    n_stages = 1
    if 'g' in fbt.keys() and not parallel:

        g_ind = dpca_obj.labels.index('g')
        g_keys = np.unique([k[g_ind] for k in plot_params.keys()])
        n_stages = len(g_keys)
        g_order = {'Cue': 0, 'PreDist': 1, 'Gating': 2, 'PostDist': 3, 'Target': 4}
        sorted_g_keys = sorted(g_keys, key=lambda x: g_order[x])
        g_keys_order = {g_key: ind for ind, g_key in enumerate(sorted_g_keys)}

        def cast_array_to_new_shape(array, order, num_stages):
            new_array = np.empty(num_stages * array.shape[0])
            new_array[:] = np.nan
            new_array[order * array.shape[0]: (order + 1) * array.shape[0]] = array
            return new_array

        for key, params in plot_params.items():
            params['line_points'] = tuple(cast_array_to_new_shape(array, g_keys_order[key[g_ind]], n_stages)
                                          for array_ind, array
                                          in enumerate(params['line_points']))


    ax2 = []
    if True:
        # init methods and params
        line_points_range = lambda line_points, ind_max: [dim_points[:ind_max] for dim_points in line_points]
        plot_line_params = lambda ax, prm: ax.plot(*line_points_range(prm['line_points'], 1), color=prm['color'],
                                                   linestyle=prm['linestyle'], linewidth=prm['linewidth'], alpha=prm['alpha'])

        def update_lines(t_i, tb_columns, plot_params, lines):
            for line, (cond, line_params) in zip(lines, plot_params.items()):
                # if cond_getter(cond) in cond_list:
                data = np.array(line_points_range(line_params['line_points'], t_i))
                line.set_data(data[:2, :])
                if len(line_params['line_points']) == 3:
                    line.set_3d_properties(data[2, :])
            return lines


        # plot
        fig = plt.figure()
        ax = plt.subplot('111') if n_dims < 3 else plt.subplot(111, projection='3d')
        # create line objects
        lines = [plot_line_params(ax, line_params)[0] for line_params in plot_params.values()]

        # Setting the axes properties
        if n_dims < 3:
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
        elif n_dims == 3:
            ax.set_xlim3d(xlim)
            ax.set_ylim3d(ylim)
            ax.set_zlim3d(zlim)
            ax.view_init(elev=elev, azim=azim)

        if legend_show:
            leghand = ax.legend(handles=legend, loc='upper right', bbox_to_anchor=(1, 1), frameon=False, prop={'size': 9})
            for legobj in leghand.legendHandles:
                legobj.set_linewidth(3.0)
        label_list = ['Time ($ms$)' if md[1] < 0 else '{0:s}-PC{1:d}'.format(md[0], md[1] + 1) for md in margin_dim_list]
        ax.set_xlabel(label_list[0])
        ax.set_ylabel(label_list[1])
        if n_dims == 3:
            ax.set_zlabel(label_list[2])
        # plt.title(' - '.join([comp[md[0]] + ' PC' + str(md[1]) for md in margin_dim_list]))

        if n_dims == 2:
            if margin_dim_list[0][1] < 0:
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)


        # set azim, elev
        # set axis labels

        # Creating the Animation object
        line_ani = FuncAnimation(fig, update_lines, len(tb_columns) * n_stages, fargs=(tb_columns, plot_params, lines),
                                 interval=50, repeat_delay=1000, blit=True)
        # plt.show()
        line_ani.save("figures/dPCA/{0:s}.gif".format(ver_str), writer=animation.PillowWriter(fps=18))


        ax2 = ax
        plt.close()


    return ax1, ax2, dpca_obj






###################
###################
###################


def dpca_pev_plot(version, num_factors, inner_show, legend_show, title_str=None):

    dpca = DemixedPrincipalComponent(DataBase([]), version)
    db, md = dpca.db, dpca.db.md

    # overwrite check
    target_filename = [dpca.get_path_base(fn, dpca.get_exec_stem()) for fn in ['fbt', 'dpca_obj', 'X_fit', 'X_tuple']]
    fbt, dpca_obj, X_fit, (X, X_trial, records, records_trial) = [md.np_loader(fn) for fn in target_filename]

    version_factor = version['factor']
    version_fr = version['fr']
    version_areas = version['area']
    version_subjects = version['subject']
    counts_thr = int(version['counts_thr'])
    ver_str = 'dPCA_{0:s}_{1:s}_{2:s}_{3:s}_{4:03d}'.format(version_factor, version_fr, version_areas, version_subjects, counts_thr)
    if title_str is None:
        title_str = 'dPCA_{0:s}_{1:s}_{2:s}\n{3:s}_{4:03d}'.format(version_factor, version_fr, version_areas, version_subjects, counts_thr)


    # plotting parameters
    color_dict, order = marginal_plotting_parameters()
    _, _, _, comp = plotting_parameters(version_factor)
    marginalizations = sorted(dpca_obj.marginalizations, key=lambda el: order.index(el))

    margin_variance = {factor_tuple: var_calc(dpca_obj, X, factor_tuple, marginalizations) for factor_tuple in
                       product(marginalizations, range(num_factors))}  # range(dpca.n_components))}
    # total_variance = {factor_tuple: sum(margin_var.values()) for factor_tuple, margin_var in margin_variance.items()}
    var_vals = np.array([list(marg.values()) for marg in margin_variance.values()])
    vals = var_vals / np.sum(var_vals)
    outer_vals = vals.sum(axis=1).flatten()
    inner_vals = vals.flatten()

    outer_colors = np.array([color_dict[m] for (m, _), margin_var in margin_variance.items()]).flatten()
    outer_marg_labels = np.array([m for (m, _), margin_var in margin_variance.items()]).flatten()
    inner_colors = np.array([color_dict[m] for _, margin_var in margin_variance.items() for m, _ in margin_var.items()]).flatten()
    size = 0.4
    radius = 1
    explode = 0



    # get margin wedge centers
    plt.ioff()
    fig, ax = plt.subplots()
    agg_vals = outer_vals.reshape(-1, num_factors).sum(axis=1).flatten()
    ax.pie(agg_vals, explode=np.array([explode for ii in range(len(agg_vals))]), radius=1,
           wedgeprops=dict(width=size, edgecolor='w'))
    margin_center_list = [ch.center for ch in ax.get_children() if type(ch) == Wedge]
    plt.close()
    plt.ion()

    # plot
    fig, ax = plt.subplots()
    pie_wedge_collection = ax.pie(outer_vals, radius=radius, colors=outer_colors,
                                                    wedgeprops=dict(width=size, edgecolor='w', linewidth=.33))
    outer_centers = np.tile(margin_center_list, num_factors).reshape(-1, 2)
    bbox_props = dict(boxstyle='square,pad=0.3', lw=0.72)
    kw = dict(arrowprops=dict(arrowstyle='-'),
              bbox=bbox_props, zorder=0, va='center')
    for ii, wedge in enumerate(pie_wedge_collection[0]):
        #coloring
        wedge.set_edgecolor('#eeeeee')  # outer_colors[ii])
        # wedge.set_alpha(1 - (ii % num_factors) * 0.2)
        wedge.set_alpha(1)
        wedge.set_center(outer_centers[ii, :])
        # text
        ang = (wedge.theta2 - wedge.theta1) / 2. + wedge.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: 'right', 1: 'left'}[int(np.sign(x))]
        connectionstyle = 'angle,angleA=0,angleB={}'.format(ang)
        kw['arrowprops'].update({'connectionstyle': connectionstyle}, color=outer_colors[ii])
        kw['bbox'].update(fc=outer_colors[ii], ec=outer_colors[ii], alpha=0.8)
        if outer_vals[ii] > .025:
            ax.annotate('{0:s}-PC{1:d} - {2:d}%'.format(outer_marg_labels[ii], (ii % num_factors) + 1, int(100 * outer_vals[ii])),
                        xy=(x, y), xytext=(1.15 * np.sign(x), 1.2 * y),
                        horizontalalignment=horizontalalignment, fontsize=8, **kw)
    if inner_show:
        pie_wedge_collection = ax.pie(inner_vals, radius=radius - size - 0.025, colors=inner_colors,
                                      wedgeprops=dict(width=0.08, edgecolor='w', linewidth=.33))
        inner_centers = np.tile(margin_center_list, num_factors * len(marginalizations)).reshape(-1, 2)
        for ii, wedge in enumerate(pie_wedge_collection[0]):
            wedge.set_edgecolor('#eeeeee')  # outer_colors[ii])
            wedge.set_center(inner_centers[ii, :])
    ax.set(aspect='equal', title='Pie plot with `ax.pie`')
    legend_elements = [Line2D([0], [0], marker='o', color='w', label=comp[m],
                              markerfacecolor=color_dict[m], markersize=12)
                       for m in marginalizations]
    if legend_show:
        ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.35, 0.9), frameon=False, prop={'size': 9})
    plt.show()
    ax.set_title(title_str)

    plt.savefig(path.join('figures', 'dPCA', '{0:s}_Variance.png'.format(ver_str)), format='png')

    return ax



################
################
################



def var_calc(dpca_obj, X_unit, factor_tuple, marginalizations):
    var = lambda Xr, D, k, tot_v: np.sum(np.dot(D[:, k], Xr) ** 2) / tot_v
    X = dpca_obj._zero_mean(X_unit)
    total_variance = np.sum((X - np.mean(X)) ** 2)
    marginal, k = factor_tuple
    D = dpca_obj.D[marginal]
    X_margin = dpca_obj._marginalize(X)
    return {m: var(X_margin[m], D, k, total_variance) for m in marginalizations}


def plotting_parameters(version_factor):

    if version_factor == 'StimulusGating':
        linestyle_column = 'StageStimSpecialized'
        color_column = 'GatingCondSpecialized'
    elif version_factor == 'StimulusGatingNoTarget':
        linestyle_column = 'StageStimSpecialized'
        color_column = 'GatingCondSpecialized'
    elif version_factor == 'StimulusGatingBool':
        linestyle_column = 'StageStimSpecialized'
    elif version_factor == 'StimulusGatingPreBool':
        linestyle_column = 'StageStimSpecialized'
        color_column = 'GatingCondSpecialized'
    elif version_factor == 'RuleStimGating':
        linestyle_column = 'RuleStimCategory'
        color_column = 'GatingCondSpecialized'
    elif version_factor == 'RuleStimGatingBool':
        linestyle_column = 'RuleStimCategory'
        color_column = 'GatingBoolSpecialized'
    elif version_factor == 'TargPostStimulusRuleStim':
        linestyle_column = 'PostStageStimSpecialized'
        color_column = 'PostRuleStimCategory'
    elif version_factor == 'GatPostStimulusRuleStim':
        linestyle_column = 'GatPostRuleStimCategory'
        color_column = 'GatPostStageStimSpecialized'
    elif version_factor == 'RuleStimGatingNull':
        linestyle_column = 'RuleStimCategory'
        color_column = 'GatingNullCondSpecialized'
    elif version_factor == 'GatedStimulus':
        linestyle_column = 'GatedStageStimSpecialized'
        color_column = 'GatedStageStimSpecialized'
    elif version_factor == 'GatedStimulusPostDistMemory':
        linestyle_column = 'GatedStimulusPostDistMemory'
        color_column = 'GatedStimulusPostDistMemory'
    elif version_factor == 'GatingPreBool':
        linestyle_column = 'GatingCondSpecialized'
        color_column = 'GatingCondSpecialized'

    alpha_column = linestyle_column

    linestyle_list, linestyle_cond_list = [], []
    color_list, color_cond_list = [], []
    alpha_list, alpha_cond_list = [], []

    if version_factor in ['GatedStimulus', 'GatedStimulusPostDistMemory']:

        linestyle_list = [(0, ()), (0, ()), (0, ()), (0, ())]
        linestyle_cond_list = ['S11', 'S12', 'S21', 'S22']

        color_list = ['#42d680', '#5e963e', '#ba58c7', '#782b43']
        color_cond_list = ['S11', 'S12', 'S21', 'S22']

        alpha_list = [1, 1, 1, 1]
        alpha_cond_list = ['S11', 'S12', 'S21', 'S22']

    elif version_factor in ['GatingPreBool']:

        linestyle_list = [(0, ()), (0, ())]
        linestyle_cond_list = ['PreDist', 'Gating']

        color_list = ['#dbb356', '#3c3a7a']
        color_cond_list = ['PreDist', 'Gating']

        alpha_list = [1, 1]
        alpha_cond_list = ['Gating', 'PreDist']

    else:
        # linestyle_list = [(0, ()), (0, (6, 2)), (0, (1, 1, 1, 1, 1, 9)), (0, (4, 7, 1, 7, 1, 7))]
        # linestyle_cond_list = ['S11', 'S12', 'S21', 'S22']

        linestyle_list = [(0, ()), (0, ()), (0, ()), (0, ())]
        linestyle_cond_list = ['S11', 'S12', 'S21', 'S22']

        alpha_list = [1, 0.8, 0.6, 0.4]
        alpha_cond_list = ['S11', 'S12', 'S21', 'S22']

        color_list = ['#42d680', '#5e963e', '#ba58c7', '#782b43',
                      '#dbb356', '#3c3a7a',
                      '#317275', '#dbb356', '#e08e4f', '#808080']
        color_cond_list = ['S11', 'S12', 'S21', 'S22',
                           'Dist', 'Gating',
                           'Cue', 'PreDist', 'PostDist', 'Target']

    comp_char = ['t', 'm', 's', 'g', 'x', 'mt', 'st', 'gt', 'sm', 'mg', 'sg', 'sgt', 'mgt', 'smt']
    comp_str = ['Time', 'Memory', 'Sensory', 'Stage', 'Miliseconds-non',
                'Memory-Time', 'Sensory-Time', 'Stage-Time',
                'Sensory-Memory', 'Memory-Stage', 'Sensory-Stage',
                'Full', 'Full', 'Full']

    linestyle = dict(zip(linestyle_cond_list, linestyle_list))
    color = dict(zip(color_cond_list, color_list))
    alpha = dict(zip(alpha_cond_list, alpha_list))
    comp = dict(zip(comp_char, comp_str))
    return (linestyle_column, linestyle), (color_column, color), (alpha_column, alpha), comp


def marginal_plotting_parameters():
    order = ['t', 's', 'st', 'm', 'mt', 'g', 'gt', 'sg', 'sm', 'mg', 'sgt', 'mgt', 'smt']
    color_dict = {'t': '#7d7d7d', 's': '#e35f6c', 'g': '#42a1eb', 'm': '#f0ed3c',
                  'st': '#9e555d', 'gt': '#457aa3', 'mt': '#b5b447',
                  'sg': '#8d4cc2', 'sm': '#c96d2a', 'mg': '#58d17c',
                  'sgt': '#694885', 'smt': '#965323', 'mgt': '#48945f'}
    return color_dict, order


def construct_margin_dim_lists():
    param = []
    for margin in ['t', 's', 'm', 'g', 'st', 'mt', 'gt', 'sm', 'sg', 'mg', 'smt', 'sgt', 'mgt']:
        param.append([(margin, 0), (margin, 1), (margin, 2)])
        param.append([(margin, 0), (margin, 1)])
        param.append([(margin, 0), (margin, 2)])
        param.append([(margin, 1), (margin, 2)])
        param.append([('x', -1), (margin, 0)])
    return param


def margin_dim_list_to_str(margin_dim_list):
    return '_'.join([m + str(i).replace('-', 'n') for m, i in margin_dim_list])



    # # SPEED OF DYNAMICS
# from statsmodels.tsa.holtwinters import SimpleExpSmoothing, Holt
# plt.figure(figsize=(16, 7))
#
# for ii in range(len(subplot_num)):
#     ax = plt.subplot(subplot_num[ii])
#     for c1 in range(C1):
#         for c2 in range(C2):
#
#                 linewidth_cc = 0.85
#                 if version_factor == 'StimulusGating' or version_factor == 'StimulusGatingBool':
#                     color_cc = color[condition_levels[1][c2]]
#                     linestyle_cc = linestyle[condition_levels[0][c1]]
#                 elif version_factor == 'RuleStimGating' or version_factor == 'RuleStimGatingBool':
#                     color_cc = color[condition_levels[1][c2]]
#                     linestyle_cc = linestyle[condition_levels[0][c1]]
#                 elif version_factor == 'PostRuleStimStimulus':
#                     color_cc = color[condition_levels[0][c1]]
#                     linestyle_cc = linestyle[condition_levels[1][c2]]
#                     linewidth_cc = 3.5 if c1 == c2 else linewidth_cc
#
#                 t_depth = 1
#                 tb_columns_speed, speed = list(zip(*[(tb_columns[t],
#                                                       np.linalg.norm(Z[comp_char[ii]][:1, c1, c2, t - t_depth] - Z[comp_char[ii]][:1, c1, c2, t]))
#                                                      for t in range(t_depth, len(tb_columns))]))
#                 model = Holt(np.asarray(speed))
#                 model._index = tb_columns_speed
#                 fit = model.fit(smoothing_level=.3, smoothing_slope=.05)
#
#                 # model = ExponentialSmoothing(np.array(speed), trend='mul', seasonal=None)
#                 # model._index = tb_columns_speed
#                 # fit = model.fit()
#
#                 ax.plot(tb_columns_speed, fit.fittedvalues,
#                             color=color_cc, linestyle=linestyle_cc, linewidth=linewidth_cc)
#     plt.title('{0:s} component'.format(comp_str[ii]))
#     plt.savefig(path.join('figures', '{0:s}_Velocity.png'.format(ver_str)), format='png')
# plt.close('all')
