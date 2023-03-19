import sys
from rec_stats import *

from seaborn import heatmap
from matplotlib import colors
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
import cmasher as cmr
from matplotlib.colors import ListedColormap


class CrossTemporal:

    def __init__(self, cv_score_obs, cv_score_shuffle, percentile=95):

        self.cv_score_obs = cv_score_obs
        self.cv_score_shuffle = cv_score_shuffle
        self.percentile = percentile

    @property
    def cv_score_thr(self):
        return np.nanpercentile(self.cv_score_shuffle, self.percentile, axis=-1)

    def get_mask(self, tail='h'):

        sc_shuffle = self.cv_score_shuffle
        sc_obs = self.cv_score_obs
        sc_thr = self.cv_score_thr
        cluster_distr, cluster_thr, cluster_obs = [], [], []

        sys.setrecursionlimit(4500)
        cluster_distr = [estimate_cluster_mass_2d(sc_shuffle[:, :, sh_i] - sc_thr, 'greater' if tail == 'h' else 'less')
                           for sh_i in range(sc_shuffle.shape[2])]
        cluster_thr = np.nanpercentile([cl_val for cl_val in cluster_distr if bool(cl_val)], 97.5)  # 97.5
        cluster_obs = estimate_cluster_mass_2d(sc_obs[:, :] - sc_thr, 'greater' if tail == 'h' else 'less', True)
        comparator = lambda x: x[0] > cluster_thr if tail == 'h' else x[0] < cluster_thr
        cluster_obs = list(filter(comparator, cluster_obs))
        sys.setrecursionlimit(3000)

        mask = coords_list_to_mask([coords for l in cluster_obs for coords in l[1]], (sc_obs.shape[0], sc_obs.shape[1]))

        return mask


def plot_cross_temporal(ax, cv_score_obs, mask):

    from seaborn import heatmap
    from matplotlib import colors
    from scipy.ndimage.filters import gaussian_filter

    t = TimebinInterval(150, 25, -200, 1000)

    def interp_t(t_ind):
        return np.interp(t_ind, t.split_to_bins_offset(), np.arange(len(t.split_to_bins_offset())))

    # def interp_inv_t(t_i):
    #     return np.interp(t_i, np.arange(len(t.split_to_bins_offset())), t.split_to_bins_offset())

    ticks_labels = [0, 385, 950]
    ticks = [interp_t(line) for line in ticks_labels]
    cmap = plt.get_cmap('inferno_r')

    def get_alpha_blend_cmap(cmap, alpha):
        cls = plt.get_cmap(cmap)(np.linspace(0, 1, 256))
        cls = (1 - alpha) + alpha * cls
        return ListedColormap(cls)

    alpha_blend = 0.6
    cmap_blend = get_alpha_blend_cmap(cmap, alpha_blend)

    heatmap(gaussian_filter(cv_score_obs, sigma=0.75), ax=ax, cmap=cmap_blend,
            norm=colors.TwoSlopeNorm(vmin=.3, vcenter=.6, vmax=.9), cbar=False, linewidths=0)
    heatmap(gaussian_filter(cv_score_obs, sigma=0.75), ax=ax, cmap=cmap,
            norm=colors.TwoSlopeNorm(vmin=.3, vcenter=.6, vmax=.9), cbar=False, mask=mask, linewidths=0)
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
