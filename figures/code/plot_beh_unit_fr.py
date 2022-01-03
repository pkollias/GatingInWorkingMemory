from rec import *
from rec_analyses import *
from rec_format import *
from versioning import *
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle



def beh_unit_fr_plot(u_ind, version_fr, version_condition, legend_show=True, title_str=None):

    # u_iloc = 0
    # version_fr = 'ConcatFactor2'
    # version_condition = 'PresentedImage'

    v_fr_params = version_fr_params(version_fr)
    t_start = v_fr_params['t_start']
    t_end = v_fr_params['t_end']
    timebin = v_fr_params['timebin']
    timestep = v_fr_params['timestep']
    timebin_interval = TimebinInterval(timebin, timestep, t_start, t_end)

    md = MetaData()
    db = md.db_base_loader(['units', 'trials', 'conditions', 'events'])
    units, trials, conditions, events = db['units'], db['trials'], db['conditions'], db['events']
    events_index = md.preproc_imports['events']['index']
    conditions.set_index(events_index, inplace=True, drop=False)
    events_conditions = pd.merge(events.reset_index(drop=True), conditions.reset_index(drop=True), on=events_index)
    events_conditions.set_index(events_index, inplace=True, drop=False)
    condition_columns, condition_series = get_condition_series(events_conditions, version_condition)

    if type(u_ind) == int:
        unit_entry = units.iloc[u_ind]
    elif type(u_ind) == tuple:
        unit_entry = units.loc[u_ind]
    unit_ind = tuple(unit_entry[md.preproc_imports['units']['index']])
    sess, channum, unitnum = unit_ind
    src_filename = md.proc_dest_path(path.join('BehavioralUnits', 'FiringRates',
                                               behunit_params_str(version_fr, timebin, timestep, t_start, t_end)),
                                     'behunit_FR_{0:s}_chan{1:03d}_unit{2:03d}.pkl'.format(sess, channum, unitnum))
    timebin_fr_dict = md.np_loader(src_filename)
    timebin_fr_dict = {key: val for key, val in timebin_fr_dict.items() if trials.loc[key[:2]].StopCondition == 1}

    pbt_mean, pbt_se = get_unit_linepoints(timebin_fr_dict, unit_ind, timebin_interval, condition_columns, condition_series)

    fig, ax = plot_unit(unit_ind, pbt_mean, pbt_se, version_condition, legend_show, title_str)
    return fig, ax




def get_unit_linepoints(timebin_fr_dict, unit_ind, timebin_interval, condition_columns, condition_series):

    # cond_val_from_entry = lambda x: tuple(x) if len(x) > 1 else x[0]
    pbt = PopulationBehavioralTimeseries(condition_columns, timebin_interval)
    data_list = []
    for event_ind in timebin_fr_dict.keys():
        entry_val = condition_series.loc[event_ind] ####### STAYED HERE
        if (type(entry_val) == float and not np.isnan(entry_val)) or type(entry_val) != float:
            cond_val = entry_val #cond_val_from_entry(entry_val)
            data_list.append([unit_ind, event_ind, cond_val, np.nan, timebin_fr_dict[event_ind]])
    pbt.add_data_rows_from_list(data_list)
    smoother = SignalSmoothing(signal.correlate, signal.windows.gaussian(13, 9))
    pbt_smooth = pbt.init_with_df(pbt.smooth_df(smoother))
    pbt_smooth.crop_timeseries(-50 - timebin_interval.timebin, 950)

    pbt_mean = pbt_smooth.average_instances(['Unit', 'Condition'])
    pbt_se = pbt_smooth.se_instances(['Unit', 'Condition'])

    return pbt_mean, pbt_se



def get_condition_series(events_conditions, version_condition):

    if version_condition == 'PresentedImage':
        condition_columns = ['StageStimExtended']
        condition_series = events_conditions[condition_columns[0]]
    elif version_condition == 'Stage':
        condition_columns = ['GatingCondExtended']
        condition_series = events_conditions[condition_columns[0]]
    elif version_condition == 'StimulusPreDistGating':
        condition_columns = ['GatingCondExtended', 'StageStimExtended']
        condition_series = events_conditions[condition_columns].apply(lambda row: tuple(row) if row['GatingCondExtended'] in ['PreDist', 'Gating'] and
                                                                                                row['StageStimExtended'] in ['S11', 'S12', 'S21', 'S22'] else np.nan, axis=1)
    elif version_condition == 'GatedStimulusGatingPostDist':
        condition_columns = ['GatingCondExtended', 'RuleStimCategory']
        condition_series = events_conditions[condition_columns].apply(lambda row: tuple(row) if row['GatingCondExtended'] in ['Gating', 'PostDist'] and
                                                                                                row['RuleStimCategory'] in ['S11', 'S12', 'S21', 'S22'] else np.nan, axis=1)
    elif version_condition == 'PresentedStimulus':
        condition_columns = ['StageStimExtended']
        condition_series = events_conditions[condition_columns[0]].replace(['C11', 'C12', 'C21', 'C22', 'S00'], np.nan)
    elif version_condition == 'GatedStimulus':
        sel_cols = ['RuleStimCategory', 'GatingCondExtended']
        condition_columns = ['RuleStimCategory']
        condition_series = events_conditions[sel_cols].apply(lambda row: row['RuleStimCategory'] if row['GatingCondExtended'] in ['Gating'] and
                                                                                                    row['RuleStimCategory'] in ['S11', 'S12', 'S21', 'S22'] else np.nan, axis=1)
    elif version_condition == 'GatedStimulusPostDist':
        sel_cols = ['RuleStimCategory', 'GatingCondExtended']
        condition_columns = ['RuleStimCategory']
        condition_series = events_conditions[sel_cols].apply(lambda row: row['RuleStimCategory'] if row['GatingCondExtended'] in ['PostDist'] and
                                                                                                    row['RuleStimCategory'] in ['S11', 'S12', 'S21', 'S22'] else np.nan, axis=1)
        condition_series

    return condition_columns, condition_series



def get_condition_plot_params(condition, version_condition):

    color_list = ['#70b349', '#518235', '#cc68d9', '#9c50a6', '#888888',
                  # '#42d680', '#5e963e', '#ba58c7', '#782b43', '#888888',
                  '#2396cc', '#2396cc', '#cf9823', '#cf9823',
                  '#dbb356', '#3c3a7a',
                  '#317275', '#dbb356', '#e08e4f', '#808080',
                  # '#e6c243', '#b39734', '#d6aa4b', '#ad8b40',
                  # '#4946a3', '#3a3885', '#4364b0', '#2b4278',
                  '#70b349', '#518235', '#cc68d9', '#9c50a6',
                  '#70b349', '#518235', '#cc68d9', '#9c50a6']
    shade_list = color_list
    linestyle_list = [(0, ()), (0, ()), (0, ()), (0, ()), (0, ()),
                      (0, ()), (0, ()), (0, ()), (0, ()),
                      (0, ()), (0, ()),
                      (0, ()), (0, ()), (0, ()), (0, ()),
                      (0, (1, 3)), (0, (1, 3)), (0, (1, 3)), (0, (1, 3)),
                      (0, (1, 3)), (0, (1, 3)), (0, (1, 3)), (0, (1, 3)),
                      (0, ()), (0, ()), (0, ()), (0, ())]
    cond_list = ['S11', 'S12', 'S21', 'S22', 'S00',
                 'C11', 'C12', 'C21', 'C22',
                 'Dist', 'Gating',
                 'Cue', 'PreDist', 'PostDist', 'Target',
                 ('PreDist', 'S11'), ('PreDist', 'S12'), ('PreDist', 'S21'), ('PreDist', 'S22'),
                 ('PostDist', 'S11'), ('PostDist', 'S12'), ('PostDist', 'S21'), ('PostDist', 'S22'),
                 ('Gating', 'S11'), ('Gating', 'S12'), ('Gating', 'S21'), ('Gating', 'S22')]
    color = dict(zip(cond_list, color_list))
    shade_color = dict(zip(cond_list, shade_list))
    linestyle = dict(zip(cond_list, linestyle_list))


    # Generalized Gating mod
    if version_condition in ['StimulusPreDistGating', 'GatedStimulusGatingPostDist']:
        clr_condition = condition[0]
        shd_clr_condition = condition[0]
        ls_condition = condition[0]
    else:
        clr_condition = condition
        shd_clr_condition = condition
        ls_condition = condition

    #
    # clr_condition = condition
    # shd_clr_condition = condition
    # ls_condition = condition

    return {'color': color[clr_condition],
            'shade_color': shade_color[shd_clr_condition],
            'linestyle': linestyle[ls_condition]}



def plot_unit(unit_ind, pbt_mean, pbt_se, version_condition, legend_show=True, title_str=None):

    # fig = plt.figure(figsize=(6.5, 5))
    fig = plt.figure(figsize=(5.55, 4.35))
    ax = plt.subplot('111')
    tb_columns = np.array(pbt_mean.timebin_interval.split_to_bins_offset())
    for row_mean, row_se in zip(pbt_mean.df.itertuples(), pbt_se.df.itertuples()):
        cond = row_mean.Condition
        sig_mean = np.array(row_mean.Timeseries)
        sig_se = np.array(row_se.Timeseries)
        plot_params = get_condition_plot_params(cond, version_condition)
        clr = plot_params['color']
        shd_clr = plot_params['shade_color']
        ls = plot_params['linestyle']
        ax.plot(tb_columns, sig_mean, color=clr, linewidth=2.5, linestyle=ls)
        ax.fill_between(tb_columns, sig_mean - sig_se, sig_mean + sig_se, alpha=0.1, edgecolor=shd_clr, facecolor=shd_clr)

    if legend_show:
        ax.legend(list(pbt_mean.df['Condition']), loc='upper right', bbox_to_anchor=(1, 1), frameon=False, prop={'size': 9})

    # ax.set_ylim(-0.15711426362251285, 12.028250268289863)

    ylims = ax.get_ylim()
    ax.add_patch(Rectangle((0, ylims[0]), 385, ylims[1] - ylims[0], color='k', alpha=0.1))
    if title_str is None:
        title_str = '{0:s} - Ch{1:03d} - Un{2:03d}'.format(*unit_ind)
    ax.set_title(title_str)
    ax.set_xlabel('event time ($ms$) relative to sample onset')
    ax.set_ylabel('firing rate ($Hz$)')

    ax.set_position([0.12, 0.1, 0.6, 0.8])
    if legend_show:
        ax.get_legend().set_bbox_to_anchor((1, .5))
        ax.get_legend()._set_loc(6)
    ax.autoscale(enable=True, axis='x', tight=True)

    t = TimebinInterval(150, 25, -50 - 150, 1000)
    # def interp_t(t_ind): return np.interp(t_ind, t.split_to_bins_offset(), np.arange(len(t.split_to_bins_offset())))
    ticks_labels = [0, 385, 950]
    ax.set_xticks(ticks_labels)
    ax.set_xticklabels(ticks_labels)

    plt.box(False)

    return fig, ax
