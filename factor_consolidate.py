import sys
from rec import *
from rec_stats import *
from factor_format import *


def main():


    # load analysis parameters
    args = sys.argv
    version_factor = args[1]
    version_fr = args[2]
    version_filter = args[3]
    version_areas = args[4]
    version_subjects = args[5]
    counts_thr = int(args[6])
    factor_method = args[7]



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
                                               factor_filter_params_str(filter_params_list_str, counts_thr, units_str), 'results'),
                                     'factor_results.pkl')
    target_filename = [md.proc_dest_path(path.join('BehavioralUnits', 'Factorization', version_factor,
                                                  behunit_params_str(version_fr, timebin, timestep, t_start, t_end),
                                                  factor_filter_params_str(filter_params_list_str, counts_thr, units_str), 'consolidate', factor_method),
                                        'factor_results_obj.pkl'),
                       md.proc_dest_path(path.join('BehavioralUnits', 'Factorization', version_factor,
                                                   behunit_params_str(version_fr, timebin, timestep, t_start, t_end),
                                                   factor_filter_params_str(filter_params_list_str, counts_thr, units_str), 'consolidate', factor_method),
                                         'decomp_obj.pkl'),
                       md.proc_dest_path(path.join('BehavioralUnits', 'Factorization', version_factor,
                                                   behunit_params_str(version_fr, timebin, timestep, t_start, t_end),
                                                   factor_filter_params_str(filter_params_list_str, counts_thr, units_str), 'consolidate', factor_method),
                                         'X_unit.pkl')]
    print(target_filename)
    if all([path.exists(fn) for fn in target_filename]):
        exit()



    factor_results = md.np_loader(src_filename)

    condition_labels = factor_results['condition_labels']
    timebin_interval = factor_results['timebin_interval']
    if factor_method == 'dPCA':
        records = factor_results['dpca_dict']['records']
        decomp_obj = factor_results['dpca_dict']['decomp_obj']
        X_factor = {margin: x for margin, x in factor_results['dpca_dict']['X_factor'].items()}
        X_unit = factor_results['dpca_dict']['X_unit']
        fbt = {margin: FactorBehavioralTimeseries(fbt_df_from_dPCA(xf, records, decomp_obj, timebin_interval), condition_labels, timebin_interval)
               for margin, xf
               in X_factor.items()}
    elif factor_method == 'PCA_mean':
        records = factor_results['pca_mean_dict']['records']
        decomp_obj = factor_results['pca_mean_dict']['decomp_obj']
        X_factor = factor_results['pca_mean_dict']['X_factor']
        X_unit = factor_results['pca_mean_dict']['X_unit']
        fbt = FactorBehavioralTimeseries(fbt_df_from_PCA(X_factor, records, decomp_obj, timebin_interval), condition_labels, timebin_interval)
    elif factor_method == 'PCA_instance':
        records = factor_results['pca_instance_dict']['records']
        decomp_obj = factor_results['pca_instance_dict']['decomp_obj']
        X_factor = factor_results['pca_instance_dict']['X_factor']
        X_unit = factor_results['pca_instance_dict']['X_unit']
        fbt = FactorBehavioralTimeseries(fbt_df_from_PCA(X_factor, records, decomp_obj, timebin_interval), condition_labels, timebin_interval)


    for f, fn in zip([fbt, decomp_obj, X_unit], target_filename):
        md.np_saver(f, fn)


main()
