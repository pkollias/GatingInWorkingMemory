####################################################################
####################################################################
####################################################################

# CROSS TEMPORAL
# ANGLE
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

angle = np.zeros((4, 9, 43, 43))

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

    fname = 'angle'
    target_filename = classifier_list[0].get_path_base(fname, stem, cross=True)
    try:
        angle[class_ravel_ind, split_ravel_ind, :, :] = md.np_loader(target_filename)
    except:
        pass


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

                # mask = coords_list_to_mask([coords for l in cluster_h_obs[class_to_ind, split_to_ind] for coords in l[1]], (43, 43))
                heatmap(angle[class_to_ind, split_to_ind], ax=ax, cmap=cmap,
                        norm=colors.TwoSlopeNorm(vmin=0, vcenter=60, vmax=90),
                        cbar=False, linewidths=0)
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
# path ='figures/Classification'
fname = '/CrossTemporal_{0:s}_{1:s}_{2:s}_Angle'.format(area, classifier.version['split'],
                                                             'PrePost' if 'PrePost' in version['balance_list']
                                                             else 'Centered')
format = 'svg'


ax.get_figure().savefig(path + fname + '.' + format, format=format)
