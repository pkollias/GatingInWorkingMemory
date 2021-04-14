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
    v_fr_params = anova_version_fr_params(version_fr)
    t_start = v_fr_params['t_start']
    t_end = v_fr_params['t_end']
    timebin = v_fr_params['timebin']
    timestep = v_fr_params['timestep']

    v_filter_params = factor_version_filter_params(version_filter)
    filter_params_list_str = v_filter_params['filter_params_list_str']

    v_units_area_list = factor_version_units_area_list(version_areas)
    area_list_str = v_units_area_list['area_list_str']

    v_units_subject_list = factor_version_units_subject_list(version_subjects)
    subject_list_str = v_units_subject_list['subject_list_str']

    units_str = '_'.join([area_list_str, subject_list_str])



    # load data and init vars, tables, and slices
    md = MetaData()

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
    smoother = SignalSmoothing(signal.correlate, signal.windows.gaussian(9, 1.25))
    if version_fr == 'ConcatFactor2':
        pbt_exec = pbt.init_with_df(pbt.smooth_df(smoother))
        pbt_exec.crop_timeseries(-50 - timebin, 1000)
    else:
        pbt_exec = pbt


    # PCA instance
    X_unit_pca_instance, records_pca_instance = pbt_exec.to_PCA_array()
    X_unit_pca_instance = StandardScaler().fit_transform(X_unit_pca_instance.transpose()).transpose()
    pca = PCA()
    X_factor_pca_instance = pca.fit_transform(X_unit_pca_instance)
    pca_instance_dict = {'X_unit': X_unit_pca_instance,
                         'X_factor': X_factor_pca_instance,
                         'decomp_obj': pca,
                         'records': records_pca_instance}

    # PCA_mean
    X_unit_pca_mean, records_pca_mean = pbt_exec.average_instances(['Unit', 'Condition']).to_PCA_array()
    X_unit_pca_mean = StandardScaler().fit_transform(X_unit_pca_mean.transpose()).transpose()
    pca = PCA()
    X_factor_pca_mean = pca.fit_transform(X_unit_pca_mean)
    pca_mean_dict = {'X_unit': X_unit_pca_mean,
                     'X_factor': X_factor_pca_mean,
                     'decomp_obj': pca,
                     'records': records_pca_mean}

    # dPCA
    X_unit_dpca_instance, records_dpca_instance = pbt_exec.to_dPCA_trial_array()
    X_unit_dpca_mean, records_dpca_mean = pbt_exec.to_dPCA_mean_array()
    mean_shape = X_unit_dpca_mean.shape
    X_unit_dpca_mean_2d = X_unit_dpca_mean.reshape((mean_shape[0], -1))
    X_unit_dpca_mean_demean = StandardScaler(with_std=False).fit_transform(X_unit_dpca_mean_2d.transpose()).transpose().reshape(mean_shape)
    dpca = dPCA.dPCA(labels=factor_dpca_labels_mapping(version_factor), regularizer='auto')
    dpca.n_trials = counts_thr
    dpca.protect = ['t']
    X_factor_dpca_mean = dpca.fit_transform(X_unit_dpca_mean_demean, X_unit_dpca_instance)
    dpca_dict = {'X_unit_instance': X_unit_dpca_instance,
                 'X_unit': X_unit_dpca_mean,
                 'X_factor': X_factor_dpca_mean,
                 'decomp_obj': dpca,
                 'records_instance': records_dpca_instance,
                 'records': records_dpca_mean}



    factor_results = {'pca_instance_dict': pca_instance_dict,
                      'pca_mean_dict': pca_mean_dict,
                      'dpca_dict': dpca_dict,
                      'condition_labels': pbt_exec.condition_labels,
                      'timebin_interval': pbt_exec.timebin_interval}

    md.np_saver(factor_results, target_filename)



main()


