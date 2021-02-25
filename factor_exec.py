import sys
from rec import *
from rec_stats import *
from factor_format import *
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from dPCA import dPCA


def main():


    # load analysis parameters
    args = sys.argv
    version_factor = args[1]
    version_fr = args[2]
    version_filter = args[3]
    version_areas = args[4]
    version_subjects = args[5]
    counts_thr = int(args[6])



    # version parameters
    v_factor_params = factor_generate_conditions(version_factor)
    condition_list = v_factor_params['condition_list']

    v_fr_params = anova_version_fr_params(version_fr)
    t_start = v_fr_params['t_start']
    t_end = v_fr_params['t_end']
    timebin = v_fr_params['timebin']
    timestep = v_fr_params['timestep']

    v_filter_params = factor_version_filter_params(version_filter)
    balance = v_filter_params['balance']
    average = v_filter_params['average']
    if not (balance or average):
        print('Data needs to be either balanced or averaged', file=sys.stderr)
        exit()
    filter_params_list_str = v_filter_params['filter_params_list_str']

    v_units_area_list = factor_version_units_area_list(version_areas)
    area_list = v_units_area_list['area_list']
    area_list_str = v_units_area_list['area_list_str']

    v_units_subject_list = factor_version_units_subject_list(version_subjects)
    subject_list = v_units_subject_list['subject_list']
    subject_list_str = v_units_subject_list['subject_list_str']

    units_str = '_'.join([area_list_str, subject_list_str])



    # load data and init vars, tables, and slices
    md = MetaData()
    db = md.db_base_loader(['units'])
    units = db['units']

    src_filename = md.proc_dest_path(path.join('BehavioralUnits', 'Factorization', version_factor,
                                               behunit_params_str(version_fr, timebin, timestep, t_start, t_end),
                                               factor_filter_params_str(filter_params_list_str, counts_thr, units_str), 'filter'),
                                     'conditions_events_filter_obj.pkl')
    target_filename = md.proc_dest_path(path.join('BehavioralUnits', 'Factorization', version_factor,
                                                  behunit_params_str(version_fr, timebin, timestep, t_start, t_end),
                                                  factor_filter_params_str(filter_params_list_str, counts_thr, units_str), 'results'),
                                        'factor_results.pkl')
    print(target_filename)
    if path.exists(target_filename):
        exit()

    pbt = md.np_loader(src_filename)

    smooth_data = np.empty(pbt.data.shape)
    smoother = SignalSmoothing(signal.correlate, signal.windows.gaussian(4, 1))
    for index in list(product(*[list(range(el)) for el in pbt.data.shape[:3]])):
        smooth_data[index] = smoother.smoothen(pbt.data[index])
    pbt.set_data(smooth_data)


    pbt_pca = pbt.base_to_PCA()
    X_pre_pca = StandardScaler().fit_transform(pbt_pca.data.transpose()).transpose()
    pca = PCA()
    X_pca = pca.fit_transform(X_pre_pca)
    pca_dict = {'X_pre_pca': X_pre_pca,
                'X_pca': X_pca,
                'pca': pca,
                'pbt_pca': PopulationBehavioralTimeseries(condition_levels=pbt_pca.condition_levels,
                                                          condition_labels=pbt_pca.condition_labels,
                                                          timebins=pbt_pca.timebins),
                'pca_base_shape': pbt_pca.base_shape}


    if average:
        dpca_dict = {}
    else:
        pbt_dpca = pbt.base_to_dPCA()
        X_pre_dpca = pbt_dpca.data
        X_pre_dpca_mean = pbt_dpca.average_instances().get_data_without_instance_dim()
        mean_shape = X_pre_dpca_mean.shape
        X_pre_dpca_mean_2d = X_pre_dpca_mean.reshape((mean_shape[0], -1))
        X_pre_dpca_demean = StandardScaler(with_std=False).fit_transform(X_pre_dpca_mean_2d.transpose()).transpose().reshape(mean_shape)
        dpca = dPCA.dPCA(labels=factor_dpca_labels_mapping(version_factor), regularizer='auto')
        dpca.protect = ['t']
        X_dpca = dpca.fit_transform(X_pre_dpca_demean, X_pre_dpca)
        dpca_dict = {'X_pre_dpca_demean': X_pre_dpca_demean,
                     'X_pre_dpca': X_pre_dpca,
                     'X_dpca': X_dpca,
                     'dpca': dpca,
                     'pbt_dpca': PopulationBehavioralTimeseries(condition_levels=pbt_dpca.condition_levels,
                                                                condition_labels=pbt_dpca.condition_labels,
                                                                timebins=pbt_dpca.timebins),
                     'dpca_base_shape': pbt_dpca.base_shape}

    factor_results = {'pca': pca_dict,
                      'dpca': dpca_dict}

    md.np_saver(factor_results, target_filename)



main()

