from rec_utils import *
from rec_stats import *
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.patches import Circle
from scipy.stats import sem
import numpy as np


def interp_t(t_ind): return np.interp(t_ind, t.split_to_bins_offset(), np.arange(len(t.split_to_bins_offset())))


def interp_inv_t(t_i): return np.interp(t_i, np.arange(len(t.split_to_bins_offset())), t.split_to_bins_offset())


t = TimebinInterval(150, 25, -50 - 150, 1000)
tb_columns = t.split_to_bins_offset()
legend_show = True

color_list = ['#eb9873', '#7eb4cc', '#6f915c']
shade_color_list = ['#ffb18c', '#9dd3eb', '#99c282']
linewidth_list = [2.5, 2.5, 2.5]
linestyle_list = [(0, ()), (0, ()), (0, ())]
line_list = ['PFC', 'Stri', 'IT']
color = dict(zip(line_list, color_list))
shade_color = dict(zip(line_list, shade_color_list))
linewidth = dict(zip(line_list, linewidth_list))
linestyle = dict(zip(line_list, linestyle_list))

for split in ['OneStimTest']:

    score_mean = {}
    score_se = {}
    legend = []
    fig = plt.figure(figsize=(5.55, 4.35))
    ax = fig.add_subplot(1, 1, 1)

    for area in ['PFC', 'Stri', 'IT']:
        print(area, split)

        args_version_list = []

        args_class = ['class=GatingPreBool']
        args_balance = ['balance=Stimulus']
        args_fr = ['fr=ConcatFactor2']
        args_counts_thr = ['counts_thr=15']
        args_area_list = ['area_list={0:s}'.format(area)]
        args_subject = ['subject=Gonzo_Oscar']
        args_area = ['area={0:s}'.format(area)]
        args_mode = ['mode=Bootstrap']
        args_mode_seed = ['mode_seed={0:d}'.format(ii) for ii in range(5)]
        args_pseudo_seed = ['pseudo_seed=0']
        args_split = ['split={0:s}'.format(split)]
        args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr,
                                                             args_area_list, args_subject, args_area, args_mode,
                                                             args_mode_seed, args_pseudo_seed, args_split)))))

        args_version_list = args_version_list
        cv_score = np.empty((43, 43, len(args_version_list)))
        for ii, args_version in enumerate(args_version_list):
            # load analysis parameters
            print(ii)
            version = job_scheduler(args_version)

            # create analysis object
            classifier = ClassificationAnalysis(DataBase([]), version)
            db, md = classifier.db, classifier.db.md
            cv_score[:, :, ii] = md.np_loader(classifier.get_path_base('cv_score', classifier.get_train_test_stem()))

        from scipy.stats import sem
        from seaborn import heatmap
        from matplotlib import colors
        import matplotlib.pyplot as plt
        from scipy.ndimage.filters import gaussian_filter

        kw = {'axis': 2}
        score_mean[area] = np.mean(cv_score, **kw)
        score_se[area] = sem(cv_score, **kw)

        #
        sig_mean = score_mean[area].diagonal()
        sig_se = score_se[area].diagonal()

        ax.plot(tb_columns, sig_mean, color=color[area], ls=linestyle[area], lw=linewidth[area], alpha=1)
        ax.fill_between(tb_columns, sig_mean - sig_se, sig_mean + sig_se, alpha=0.2,
                        edgecolor=shade_color[area],
                        facecolor=shade_color[area])

        legend_str = '{0:s}'.format(area)
        # legend_str = '{1:s} (N={2:d}, {3:d}%)'.format(area, subdir_str, significant_n,
        #                                                     round(100 * (float(significant_n) / float(area_n))))
        legend_element = Circle((0, 0), radius=1, edgecolor='w', facecolor=color[area], label=legend_str)
        legend.append(legend_element)

    if legend_show:
        ax.legend(handles=legend, loc='upper right', bbox_to_anchor=(1, 1), frameon=False, prop={'size': 9})

        # ylims = ax.get_ylim()
    ylims = (0.4, 1)
    plt.ylim(ylims[0], ylims[1])
    xlims = ax.get_xlim()
    ax.add_patch(Rectangle((0, ylims[0]), 385, ylims[1] - ylims[0], color='k', alpha=0.1))
    ax.hlines(.5, tb_columns[0], tb_columns[-1], '#999999', ':')
    # ax.set_title(title_str)

    ax.set_xticks(list(np.linspace(0, 900, 4)))
    ax.set_xlabel('event time (ms)')
    ax.set_ylabel('Mean CV classification accuracy\n(averaged for cue groups)')
    ax.set_position([0.18, 0.1, 0.54, 0.8])
    if legend_show:
        ax.get_legend().set_bbox_to_anchor((1, .5))
        ax.get_legend()._set_loc(6)
    ax.autoscale(enable=True, axis='x', tight=True)
    # plt.tight_layout()
    plt.box(False)
