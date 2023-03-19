####################################################################
####################################################################
####################################################################

# CROSS TEMPORAL

from rec_utils import *
from rec_stats import *

from seaborn import heatmap
from matplotlib import colors
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
import cmasher as cmr
from matplotlib.colors import ListedColormap


args_version_list = []

area_list = 'PFC'
area = 'PFC'
args_class_list = ['class_list=Stimulus_GatedStimulus']
args_balance = ['balance_list=StageGatingPrePostMemory_StageGatingPrePostSensory']
args_fr = ['fr=ConcatFactor2']
args_counts_thr = ['counts_thr=12']
args_area_list = ['area_list={0:s}'.format(area_list)]
args_subject = ['subject=Gonzo_Oscar']
args_area = ['area={0:s}'.format(area)]
args_mode = ['mode=Normal']
args_mode_seed = ['mode_seed=0']  # 100
args_pseudo_seed = ['pseudo_seed=0']
args_split = ['split=StratifiedBalanceSplit_StimHalf']
args_classifier_ind = ['classifier_ind={0:d}_{1:d}'.format(*t) for t in list(product(range(2), range(2)))]
args_split_split_ind = ['split_split_ind={0:d}_{1:d}'.format(*t) for t in list(product(range(3), range(3)))]
args_shuffle = ['shuffle=0']
args_pbt_src = ['pbt_src=pca']
args_version_list.extend(list(map(list, list(product(args_class_list, args_balance, args_fr, args_counts_thr,
                                                     args_area_list, args_subject, args_area,
                                                     args_mode, args_mode_seed, args_pseudo_seed, args_split,
                                                     args_classifier_ind, args_split_split_ind, args_shuffle,
                                                     args_pbt_src)))))

args_version = args_version_list[0]
version = job_scheduler(args_version)
classifier = ClassificationAnalysis(DataBase([]), version)
db, md = classifier.db, classifier.db.md

cv_score_obs = np.zeros((4, 9, 43, 43))

for ii, args_version in enumerate(args_version_list):
    # load analysis parameters
    print(ii)
    version = job_scheduler(args_version)

    # create analysis object
    classifier = ClassificationAnalysis(DataBase([]), version)
    db, md = classifier.db, classifier.db.md
    classifier_ind = list(map(int, version['classifier_ind'].split('_')))
    class_ravel_ind = np.ravel_multi_index(classifier_ind, (2, 2))
    split_split_ind = list(map(int, version['split_split_ind'].split('_')))
    split_ravel_ind = np.ravel_multi_index(split_split_ind, (3, 3))

    class_list = version['class_list'].split('_')
    balance_list = version['balance_list'].split('_')

    def version_modifier(class_i, balance_i):
        version_mod = version.copy()
        version_mod['class'] = class_i
        version_mod['balance'] = balance_i
        return version_mod

    # create analysis object
    classifier_list = [ClassificationAnalysis(DataBase([]), version_modifier(class_i, balance_i)) for class_i, balance_i
                       in zip(class_list, balance_list)]
    score = {}
    # get training classifier
    class_train = classifier_list[classifier_ind[0]]
    db_train, md_train = class_train.db, class_train.db.md
    train_dict = md_train.np_loader(class_train.get_path_base('train_test', class_train.get_train_test_stem()))
    # get split of splits
    if class_train.version['split'] == 'StratifiedBalanceSplit_StimHalf':
        train_key = list(np.unique(list(map(itemgetter(0), train_dict.keys()))))[split_split_ind[0]]
    else:
        train_key = list(train_dict.keys())[split_split_ind[0]]

    # get testing classifier
    class_test = classifier_list[classifier_ind[1]]
    db_test, md_test = class_test.db, class_test.db.md
    test_dict = md_test.np_loader(class_test.get_path_base('train_test', class_test.get_train_test_stem()))
    # get split of splits
    if class_test.version['split'] == 'StratifiedBalanceSplit_StimHalf':
        test_key = list(np.unique(list(map(itemgetter(0), test_dict.keys()))))[split_split_ind[1]]
    else:
        test_key = list(test_dict.keys())[split_split_ind[1]]

    stem = (*classifier_list[0].get_train_test_stem(),
            'intermediate',
            '_'.join([class_train.version['class'], train_key, class_test.version['class'], test_key]),
            'pca' if version['pbt_src'] == 'pca' else '')

    fname = 'cv_score' if not int(version['shuffle']) else 'cv_score_{0:04d}'.format(int(version['shuffle']))
    target_filename = classifier_list[0].get_path_base(fname, stem, cross=True)
    try:
        cv_score_obs[class_ravel_ind, split_ravel_ind, :, :] = md.np_loader(target_filename)
    except:
        pass

# load shuffle data
stem = (*classifier.get_train_test_stem(), 'consolidate')
fname = 'cv_score_full'
target_filename = classifier.get_path_base(fname, stem, cross=True)
cv_score = md.np_loader(target_filename)
cv_score = cv_score[:, :, :, :, 1:]

# h and l percentiles
cv_h_thr = np.nanpercentile(cv_score, 97.5, -1)
cv_l_thr = np.nanpercentile(cv_score, 2.5, -1)

# estimate cluster mass distribution
sys.setrecursionlimit(4500)

cluster_h_distr = np.empty((cv_score.shape[:2]), dtype=list)
cluster_l_distr = np.empty((cv_score.shape[:2]), dtype=list)
cluster_h_thr = np.empty((cv_score.shape[:2]), dtype=float)
cluster_l_thr = np.empty((cv_score.shape[:2]), dtype=list)
cluster_h_obs = np.empty((cv_score.shape[:2]), dtype=list)
cluster_l_obs = np.empty((cv_score.shape[:2]), dtype=list)

for class_ind in range(cv_score.shape[0]):
    for split_ind in range(cv_score.shape[1]):
        print(class_ind, split_ind)
        inds = (class_ind, split_ind)
        cluster_h_distr[inds[0], inds[1]] = [estimate_cluster_mass_2d(cv_score[inds[0], inds[1], :, :, sh_i] - cv_h_thr[inds[0], inds[1]], 'greater')
                           for sh_i in range(cv_score.shape[4])]
        cluster_h_thr[inds[0], inds[1]] = np.nanpercentile([cl_val for cl_val in cluster_h_distr[inds[0], inds[1]] if bool(cl_val)], 97.5)
        cluster_h_obs[inds[0], inds[1]] = estimate_cluster_mass_2d(cv_score_obs[inds[0], inds[1], :, :] - cv_h_thr[inds[0], inds[1]], 'greater', True)
        cluster_h_obs[inds[0], inds[1]] = list(filter(lambda x: x[0] > cluster_h_thr[inds[0], inds[1]], cluster_h_obs[inds[0], inds[1]]))

        cluster_l_distr[inds[0], inds[1]] = [estimate_cluster_mass_2d(cv_score[inds[0], inds[1], :, :, sh_i] - cv_l_thr[inds[0], inds[1]], 'less')
                           for sh_i in range(cv_score.shape[4])]
        cluster_l_thr[inds[0], inds[1]] = np.nanpercentile([cl_val for cl_val in cluster_l_distr[inds[0], inds[1]] if bool(cl_val)], 2.5)
        cluster_l_obs[inds[0], inds[1]] = estimate_cluster_mass_2d(cv_score_obs[inds[0], inds[1], :, :] - cv_l_thr[inds[0], inds[1]], 'less', True)
        cluster_l_obs[inds[0], inds[1]] = list(filter(lambda x: x[0] < cluster_l_thr[inds[0], inds[1]], cluster_l_obs[inds[0], inds[1]]))

sys.setrecursionlimit(3000)


# plotting

axsz = 40
gszn = 3
gszj = 20
span = 40
spans = np.array([0, axsz, gszn, axsz, gszn, axsz, gszj, axsz, gszn, axsz, gszn, axsz])

t = TimebinInterval(150, 25, -50 - 150, 1000)

def interp_t(t_ind): return np.interp(t_ind, t.split_to_bins_offset(), np.arange(len(t.split_to_bins_offset())))
def interp_inv_t(t_i): return np.interp(t_i, np.arange(len(t.split_to_bins_offset())), t.split_to_bins_offset())

ticks_labels = [0, 385, 950]
ticks = [interp_t(line) for line in ticks_labels]

def pos_to_gridpos(ii): return int(np.sum(spans[:2 * ii + 1]))

cmap = plt.get_cmap('inferno_r')
# cmap, norm = colors.from_levels_and_colors(np.linspace(1, 0, 8),
#                               ['#4d2c1a', '#6d030c', '#d0141d', '#fa9e7a', '#ffc664', '#a6bac8', '#ced2d0', '#271629'],
#                               extend='max')
def get_alpha_blend_cmap(cmap, alpha):
    cls = plt.get_cmap(cmap)(np.linspace(0,1,256))
    cls = (1-alpha) + alpha*cls
    return ListedColormap(cls)
alpha_blend = 0.6
cmap_blend = get_alpha_blend_cmap(cmap, alpha_blend)

# create figure
# fig = plt.figure(figsize=(10, 10))
fig = plt.figure(figsize=(5.5, 5.5))
axs = [plt.subplot2grid((np.sum(spans), np.sum(spans)), (pos_to_gridpos(ax_i), pos_to_gridpos(ax_j)), rowspan=span,
                        colspan=span) for ax_i, ax_j in product(range(6), range(6))]
class_class_split_split_to_ind = np.arange(6 * 6).reshape((2, 3, 2, 3)).transpose((0, 2, 1, 3))
# virtually iterates [sensory, memory] x [sensory, memory] x [predist, gating, postdist] x [predist, gating, postdist]
for class_train_ind in range(2):
    for class_test_ind in range(2):
        for split_split_train_ind in range(3):
            for split_split_test_ind in range(3):

                # convert above to grid quadrants (top left 0,0, bottom right 1,1)
                class_y_ind_to_ax_ind = [1, 0]
                class_train_ax_ind = class_y_ind_to_ax_ind[class_train_ind]
                class_x_ind_to_ax_ind = [0, 1]
                class_test_ax_ind = class_x_ind_to_ax_ind[class_test_ind]

                # convert above to inner quadrant y axis flip
                split_split_x_ind_to_split_split_x_ax_ind = [2, 1, 0]
                ssxi_ssxai = split_split_x_ind_to_split_split_x_ax_ind
                split_split_train_ax_ind = ssxi_ssxai[split_split_train_ind]
                split_split_test_ax_ind = split_split_test_ind

                ax_ind = class_class_split_split_to_ind[
                    class_train_ax_ind, class_test_ax_ind, split_split_train_ax_ind, split_split_test_ax_ind]
                ax = axs[ax_ind]

                split_to_ordered_split_ind = [2, 0, 1]
                split_train_ind = split_to_ordered_split_ind[split_split_train_ind]
                split_test_ind = split_to_ordered_split_ind[split_split_test_ind]
                class_to_ind = np.ravel_multi_index((class_train_ind, class_test_ind), (2, 2))
                split_to_ind = np.ravel_multi_index((split_train_ind, split_test_ind), (3, 3))

                mask = coords_list_to_mask([coords for l in cluster_h_obs[class_to_ind, split_to_ind] for coords in l[1]], (43, 43))
                heatmap(gaussian_filter(cv_score_obs[class_to_ind, split_to_ind], sigma=0.75), ax=ax, cmap=cmap_blend,
                        norm=colors.TwoSlopeNorm(vmin=.3, vcenter=.6 if 'Half' in version['split'] else 0.33, vmax=.9),
                        cbar=False, linewidths=0)
                heatmap(gaussian_filter(cv_score_obs[class_to_ind, split_to_ind], sigma=0.75), ax=ax, cmap=cmap,
                        norm=colors.TwoSlopeNorm(vmin=.3, vcenter=.6 if 'Half' in version['split'] else 0.33, vmax=.9),
                        cbar=False, mask=mask, linewidths=0)
                # 'RdYlGn' 'PRGn' 'PiYG'
                ax.invert_yaxis()
                for hline in ticks: ax.hlines(hline, *ax.get_xlim(), color='w')
                for vline in ticks: ax.vlines(vline, *ax.get_ylim(), color='w')
                # ax.set_xlabel('Testing Classification Time (ms) from Stim Onset')
                if ax_ind == 30:
                    ax.set_xticks(ticks)
                    ax.set_xticklabels(ticks_labels)
                    # ax.set_ylabel('Training Classification Time (ms) from Stim Onset')
                    ax.set_yticks(ticks)
                    ax.set_yticklabels(ticks_labels)
                else:
                    ax.set_xticks([])
                    ax.set_yticks([])
                ax.axis('square')
                # ax.set_title(area)
fig.suptitle('{0:s}\nTrainig Set Testing Set PreDist Gating PostDist Sensory Memory'.format(area))

path ='/Users/kollias/MEGA/research/gating/Thesis/Figures/Thesis_Figures/Classification'
# path ='figures/Classification'
fname = '/CrossTemporal_{0:s}_{1:s}_{2:s}_Classifier'.format(area, classifier.version['split'],
                                                             'PrePost' if 'PrePost' in version['balance_list']
                                                             else 'Centered')
format = 'svg'


ax.get_figure().savefig(path + fname + '.' + format, format=format)


####################################################################
####################################################################
####################################################################

# intersect valid_units

from rec_analyses import *
def args_from_parse_func(parse_version):
    args_version_list = []
    for class_i, balance in [('Stimulus', 'StageGatingPrePostMemory'), ('Stimulus', 'StageGatingCenteredMemory'),
                             ('GatedStimulus', 'StageGatingPrePostSensory'), ('GatedStimulus', 'StageGatingCenteredSensory')]:
        args_class = ['class={0:s}'.format(class_i)]
        args_balance = ['balance={0:s}'.format(balance)]
        args_fr = ['fr=ConcatFactor2']
        args_counts_thr = ['counts_thr=15']
        args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr)))))
    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))
    return args_version_from_job

from rec_utils import *
args_version_loop_list = ((['job_id=0', 'overwrite=True'], ['job_id=2', 'overwrite=True']),
                          (['job_id=1', 'overwrite=True'], ['job_id=3', 'overwrite=True']))
for args_version_list in args_version_loop_list:
    for counts_thr in ['12']:
        units_events_list, filename_list = [], []
        for args_version in args_version_list:
            version = job_scheduler(args_version, args_from_parse_func)
            version['counts_thr'] = counts_thr
            classifier = ClassificationAnalysis(DataBase(['trials', 'units', 'events', 'conditions']), version)
            db, md = classifier.db, classifier.db.md
            target_filename = classifier.get_path_base('valid_units', classifier.get_wrangle_stem())
            valid_units_events = md.np_loader(target_filename)
            units_events_list.append(valid_units_events)
            filename_list.append(target_filename)
        conjoint_units_events(units_events_list, filename_list)


####################################################################
####################################################################
####################################################################

# EXAMPLES 1-D CLASSIFIERS

from rec_utils import *
from rec import *
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.patches import Circle
from scipy.stats import sem
import numpy as np

# examples 1-D classifiers

def interp_t(t_ind): return np.interp(t_ind, t.split_to_bins_offset(), np.arange(len(t.split_to_bins_offset())))
def interp_inv_t(t_i): return np.interp(t_i, np.arange(len(t.split_to_bins_offset())), t.split_to_bins_offset())

score_mean_pfc, score_sem_pfc, score_mean_stri, score_sem_stri, score_mean_it, score_sem_it, _ = MetaData().np_loader(
    'files/score.pkl')
# score_mean = cv_score.mean(axis=2)
# score_sem = sem(cv_score, axis=2)

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


legend = []
fig = plt.figure(figsize=(5.55, 4.35))
ax = fig.add_subplot(1, 1, 1)
for score_mean, score_se, area in [(score_mean_pfc, score_sem_pfc, 'PFC'), (score_mean_stri, score_sem_stri, 'Stri'),
                                   (score_mean_it, score_sem_it, 'IT')]:

    # first index: class combo [0..4] Stimulus,GatedStimulus x Stimulus,GatedStimulus
    # second index: split combo [0..9] Gating,PostDist,PreDist x Gating,PostDist,PreDist
    cc_ind = 0
    ss_ind = 0
    #
    sig_mean = score_mean[cc_ind, ss_ind, :, :].diagonal()
    sig_se = score_se[cc_ind, ss_ind, :, :].diagonal()

    # t_i = 200
    # sig_mean = score_mean[cc_ind, ss_ind, int(interp_t(t_i)), :]
    # sig_se = score_se[cc_ind, ss_ind, int(interp_t(t_i)), :]

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

# ax.set_xticks(list(np.linspace(0, 900, 4)))
ax.set_xticks([0, 385, 950])
ax.set_xlabel('event time (ms)')
ax.set_ylabel('Mean CV classification accuracy\n(averaged for cue groups)')
ax.set_position([0.18, 0.1, 0.54, 0.8])
if legend_show:
    ax.get_legend().set_bbox_to_anchor((1, .5))
    ax.get_legend()._set_loc(6)
ax.autoscale(enable=True, axis='x', tight=True)
# plt.tight_layout()
plt.box(False)


####################################################################
####################################################################
####################################################################

# GENERALIZED GATING CLASSIFIER


from rec_utils import *
from rec_stats import *

from seaborn import heatmap
from matplotlib import colors
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
import cmasher as cmr
from matplotlib.colors import ListedColormap
import re

splits_list = ['WithinGroupTransfer']
area_list = ['PFC', 'Stri', 'IT']

fig, axs = plt.subplots(len(splits_list), len(area_list))
for ii, split in enumerate(splits_list):
    for jj, area_list in enumerate(area_list):
        area = area_list
        args_class = ['class=GatingPreBoolGeneralized']
        args_balance = ['balance=Stimulus']
        args_fr = ['fr=ConcatFactor2']
        args_counts_thr = ['counts_thr=12']
        args_area_list = ['area_list={0:s}'.format(area_list)]
        args_subject = ['subject=Gonzo_Oscar']
        args_area = ['area={0:s}'.format(area)]
        args_mode = ['mode=Normal']
        args_mode_seed = ['mode_seed=0']
        args_pseudo_seed = ['pseudo_seed=0']
        args_split = ['split={0:s}'.format(split)]  # ['AcrossGroupTransfer', 'OneStimTrain', StratifiedStim']
        args_shuffle = ['shuffle=0']
        args_version = list(product(args_class, args_balance, args_fr, args_counts_thr,
                                    args_area_list, args_subject, args_area,
                                    args_mode, args_mode_seed, args_pseudo_seed, args_split,
                                    args_shuffle))[0]

        version = job_scheduler(args_version)
        classifier = ClassificationAnalysis(DataBase([]), version)
        db, md = classifier.db, classifier.db.md

        cv_score_obs = np.zeros((43, 43))

        # load analysis parameters
        version = job_scheduler(args_version)

        # create analysis object
        classifier = ClassificationAnalysis(DataBase([]), version)
        db, md = classifier.db, classifier.db.md

        # create analysis object
        # get training classifier

        stem = classifier.get_train_test_stem()
        fname = 'cv_score' if not int(version['shuffle']) else 'cv_score_{0:04d}'.format(int(version['shuffle']))
        target_filename = classifier.get_path_base(fname, stem, cross=False)
        try:
            cv_score_obs[:, :] = md.np_loader(target_filename)
        except:
            pass

        # load shuffle data
        stem = classifier.get_train_test_stem()
        fname = 'cv_score_full'
        target_filename = classifier.get_path_base(fname, stem, cross=False)
        cv_score = md.np_loader(target_filename)

        # h and l percentiles
        cv_h_thr = np.nanpercentile(cv_score, 95, -1)
        cv_l_thr = np.nanpercentile(cv_score, 5, -1)

        # estimate cluster mass distribution
        sys.setrecursionlimit(4500)

        cluster_h_distr = [estimate_cluster_mass_2d(cv_score[:, :, sh_i] - cv_h_thr, 'greater')
                           for sh_i in range(cv_score.shape[2])]
        cluster_h_thr = np.nanpercentile([cl_val for cl_val in cluster_h_distr if bool(cl_val)], 97.5) ### PK IS THAT NAN RIGHT?
        cluster_h_obs = estimate_cluster_mass_2d(cv_score_obs[:, :] - cv_h_thr, 'greater', True)
        cluster_h_obs = list(filter(lambda x: x[0] > cluster_h_thr, cluster_h_obs))

        cluster_l_distr = [estimate_cluster_mass_2d(cv_score[:, :, sh_i] - cv_l_thr, 'less')
                           for sh_i in range(cv_score.shape[2])]
        cluster_l_thr = np.nanpercentile([cl_val for cl_val in cluster_l_distr if bool(cl_val)], 2.5) ### PK IS THAT NAN RIGHT?
        cluster_l_obs = estimate_cluster_mass_2d(cv_score_obs[:, :] - cv_l_thr, 'less', True)
        cluster_l_obs = list(filter(lambda x: x[0] < cluster_l_thr, cluster_l_obs))

        sys.setrecursionlimit(3000)

        ####
        ####
        ####

        # plotting

        t = TimebinInterval(150, 25, -50 - 150, 1000)


        def interp_t(t_ind):
            return np.interp(t_ind, t.split_to_bins_offset(), np.arange(len(t.split_to_bins_offset())))


        def interp_inv_t(t_i):
            return np.interp(t_i, np.arange(len(t.split_to_bins_offset())), t.split_to_bins_offset())


        ticks_labels = [0, 385, 950]
        ticks = [interp_t(line) for line in ticks_labels]


        def pos_to_gridpos(ii):
            return int(np.sum(spans[:2 * ii + 1]))


        cmap = cmr.pride


        # cmap, norm = colors.from_levels_and_colors(np.linspace(1, 0, 8),
        #                               ['#4d2c1a', '#6d030c', '#d0141d', '#fa9e7a', '#ffc664', '#a6bac8', '#ced2d0', '#271629'],
        #                               extend='max')
        def get_alpha_blend_cmap(cmap, alpha):
            cls = plt.get_cmap(cmap)(np.linspace(0, 1, 256))
            cls = (1 - alpha) + alpha * cls
            return ListedColormap(cls)


        alpha_blend = 0.5
        cmap_blend = get_alpha_blend_cmap(cmap, alpha_blend)

        # create figure
        # fig = plt.figure(figsize=(10, 10))
        # fig, ax = plt.subplots(figsize=(3.75, 3.75))
        try:
            ax = axs[ii][jj]
        except:
            ax = axs[jj]
        # virtually iterates [sensory, memory] x [sensory, memory] x [predist, gating, postdist] x [predist, gating, postdist]

        mask = coords_list_to_mask([coords for l in cluster_h_obs for coords in l[1]], (43, 43))
        heatmap(gaussian_filter(cv_score_obs, sigma=0.75), ax=ax, cmap=cmap_blend,
                        norm=colors.TwoSlopeNorm(vmin=.3, vcenter=.6, vmax=.9),
                        cbar=False, linewidths=0)
        heatmap(gaussian_filter(cv_score_obs, sigma=0.75), ax=ax, cmap=cmap,
                        norm=colors.TwoSlopeNorm(vmin=.3, vcenter=.6, vmax=.9),
                        cbar=False, mask=mask, linewidths=0)
        # 'RdYlGn' 'PRGn' 'PiYG'
        ax.invert_yaxis()
        for hline in ticks: ax.hlines(hline, *ax.get_xlim(), color='w')
        for vline in ticks: ax.vlines(vline, *ax.get_ylim(), color='w')
        # ax.set_xlabel('Testing Classification Time (ms) from Stim Onset')
        ax.set_xticks(ticks)
        ax.set_xticklabels(ticks_labels)
        # ax.set_ylabel('Training Classification Time (ms) from Stim Onset')
        ax.set_yticks(ticks)
        ax.set_yticklabels(ticks_labels)
        ax.axis('square')
        ax.set_title(area)
        if jj == 0:
            ax.set_ylabel('{0:s}'.format(re.sub(r"(\w)([A-Z])", r"\1 \2", split)))



# GENERALIZED GATING 1-D

from rec_utils import *
from rec_stats import *
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.patches import Circle
from scipy.stats import sem
import numpy as np
from scipy.ndimage.filters import gaussian_filter


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
        args_mode_seed = ['mode_seed={0:d}'.format(ii) for ii in range(1)]
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
        sig_mean = gaussian_filter(score_mean[area], 2).diagonal()
        # sig_se = gaussian_filter(score_se[area], 2).diagonal()

        ax.plot(tb_columns, sig_mean, color=color[area], ls=linestyle[area], lw=linewidth[area], alpha=1)
        # ax.fill_between(tb_columns, sig_mean - sig_se, sig_mean + sig_se, alpha=0.2,
        #                 edgecolor=shade_color[area],
        #                 facecolor=shade_color[area])

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
