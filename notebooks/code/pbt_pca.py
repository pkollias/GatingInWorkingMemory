from rec_analyses import *
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline


#########################
''' class definitions '''
#########################

class PopulationPCs:

    def __init__(self, classes=(('Stimulus_GatedStimulus', 'StageGatingPrePostMemory_StageGatingPrePostSensory'),
                                ('Stimulus', 'StageGatingPrePostMemory'),
                                ('GatedStimulus', 'StageGatingPrePostSensory')), area='PFC'):

        def version_modifier(class_i, balance_i):
            version_mod = version.copy()
            version_mod['class'] = class_i
            version_mod['balance'] = balance_i
            return version_mod
        self.pbt = {}

        version = dict((('counts_thr', 12), ('area_list', area), ('area', area), ('subject', 'Gonzo_Oscar'),
                        ('fr', 'ConcatFactor2'), ('mode', 'Normal'), ('mode_seed', '0'), ('pca_mode', 'mean'),
                        ('class_list', np.nan), ('balance_list', np.nan)))
        for class_list, balance_list in classes:
            version['class_list'] = class_list
            version['balance_list'] = balance_list

            class_list = version['class_list'].split('_')
            balance_list = version['balance_list'].split('_')
            classifier_list = [ClassificationAnalysis(DataBase([]), version_modifier(class_i, balance_i))
                               for class_i, balance_i
                               in zip(class_list, balance_list)]
            pbt_list = [MetaData().np_loader(class_i.get_path_base('pbt', class_i.get_assemble_stem())) for class_i in
                        classifier_list]
            self.pbt[version['class_list']] = pbt_list[0].init_with_df(pd.concat([pbt_i.df for pbt_i in pbt_list]))

    def crop_timeseries(self, t_start, t_end, class_i):
        pbt_i = self.pbt[class_i]
        pbt_i.crop_timeseries(t_start, t_end)
        return self

    def average_timeseries(self, class_i):
        pbt_i = self.pbt[class_i]
        pbt_i.df['Timeseries'] = pbt_i.df['Timeseries'].apply(lambda signal: [np.mean(signal)])
        pbt_i.timebin_interval = TimebinInterval(950 - 600 + 150, 25, 600 - 150, 950)
        return self

    def slice_rename_by_cond(self, class_i, stage_list):
        pbt_i = self.pbt[class_i]
        pbt_i.df = pbt_i.df.loc[pbt_i.df['Condition'].apply(lambda cond: cond[1] in stage_list)]
        rename_lambda = lambda cond: (cond[0], '*', '*')
        pbt_i.df['Condition'] = pbt_i.df['Condition'].apply(rename_lambda)
        return self

    def average_instances(self, class_i):
        self.pbt[class_i] = self.pbt[class_i].average_instances(['Unit', 'Unit_Code', 'Condition'])
        return self

    def pipeline_fit(self, class_i):
        pbt_i = self.pbt[class_i]
        X, records = pbt_i.to_PCA_array()
        pca, scaler = PCA(), StandardScaler()
        pca.fit(scaler.fit_transform(X).T).T
        pipeline = Pipeline([('scaler', scaler), ('pca', pca)])
        return pipeline

    def pipeline_transform(self, class_i, pipeline):
        pbt_i = self.pbt[class_i]
        X, records = pbt_i.to_PCA_array()
        X_factor = pipeline.transform(X)
        return X_factor, records


def plot_trajectories(space, activity, pc_activity, timebin_interval=None, line_params=None, ax=None):

    if timebin_interval is None:
        timebin_interval = TimebinInterval(-150, 25, 150, 950)

    # space = (('Stimulus', 'PreDist', 0), ('Stimulus', 'PreDist', 1), ('Stimulus', 'PreDist', 2))
    # activity = [('Stimulus', 'Gating', 0), ('Stimulus', 'Gating', 1),
    #                ('Stimulus', 'PreDist', 2), ('Stimulus', 'PreDist', 3)]
    color_dict = {0: '#70b349', 1: '#518235', 2: '#cc68d9', 3: '#9c50a6'}
    linestyle_dict = {'PreDist': ':', 'Gating': '-', 'PostDist': '--'}
    if line_params is None:
        line_params = [{'linewidth': 2,
                        'color': color_dict.get(act[2], '#888888'),
                        'linestyle': linestyle_dict.get(act[1], '-')}
                       for act in activity]

    # plot specs
    if True:
        plot_method = lambda ax, n_dims: ax.plot if n_dims < 3 else ax.plot3D
        scatter_method = lambda ax, n_dims: ax.scatter if n_dims < 3 else ax.scatter3D
        stim_on_t = timebin_interval.split_to_bins_offset().index(0)
        stim_enc_t = timebin_interval.split_to_bins_offset().index(200)
        stim_off_t = timebin_interval.split_to_bins_offset().index(400)
        delay_off_t = timebin_interval.split_to_bins_offset().index(950)

    # plot params
    n_dims = len(space)

    # plot
    if ax is None:
        fig = plt.figure()
        ax = plt.axes(projection='3d' if n_dims == 3 else None)
    else:
        fig = ax.get_figure()
    ax.set_xlabel('{0:s}-{1:s} PC{2:d}'.format(space[0][0][0], ''.join([st[:2] for st in space[0][1].split('_')]), space[0][2] + 1))
    ax.set_ylabel('{0:s}-{1:s} PC{2:d}'.format(space[1][0][0], ''.join([st[:2] for st in space[0][1].split('_')]), space[1][2] + 1))
    if n_dims == 3: ('{0:s}-{1:s} PC{2:d}'.format(space[2][0][0], ''.join([st[:2] for st in space[0][1].split('_')]), space[2][2] + 1))
    # ax.set_title(f"{version['class_list']} - {version['stage_list']}")

    # {(pca_space, activity)} = [factor, condition, timebin]
    for line_ind, line_i in enumerate(activity):
        data_sig = [pc_activity[((pc_i[0], pc_i[1]),
                                 (line_i[0], line_i[1]))][pc_i[2], line_i[2], :]
                    for pc_ind, pc_i in enumerate(space)]
        plot_method(ax, n_dims)(*data_sig, **line_params[line_ind])
        scatter_method(ax, n_dims)(*[d_i[stim_on_t] for d_i in data_sig], color=line_params[line_ind]['color'], s=100, marker='*')
        scatter_method(ax, n_dims)(*[d_i[stim_enc_t] for d_i in data_sig], color=line_params[line_ind]['color'], s=100, marker='+')
        scatter_method(ax, n_dims)(*[d_i[stim_off_t] for d_i in data_sig], color=line_params[line_ind]['color'], s=100, marker='^')
        scatter_method(ax, n_dims)(*[d_i[delay_off_t] for d_i in data_sig], color=line_params[line_ind]['color'], s=30, marker='o')
        # scatter_method(ax, n_dims)(*data_pc[ii, :], color=color[ii], s=150, alpha=0.2)
    ax.axis('square')
    # plt.show()

    return fig, ax
