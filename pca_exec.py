import sys
from rec import *
from pca_format import *
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from dPCA import dPCA


def main():


    # load analysis parameters
    args = sys.argv
    version_pca = args[1]
    version_fr = args[2]
    version_filter = args[3]
    version_units = args[4]
    counts_thr = int(args[5])



    # version parameters
    v_pca_params = pca_generate_conditions(version_pca)
    condition_list = v_pca_params['condition_list']

    v_fr_params = anova_version_fr_params(version_fr)
    t_start = v_fr_params['t_start']
    t_end = v_fr_params['t_end']
    timebin = v_fr_params['timebin']
    timestep = v_fr_params['timestep']

    v_filter_params = pca_version_filter_params(version_filter)
    balance = v_filter_params['balance']
    average = v_filter_params['average']
    if not (balance or average):
        print('Data needs to be either balanced or averaged', file=sys.stderr)
        exit()
    filter_params_list_str = v_filter_params['filter_params_list_str']

    v_units_area_list = pca_version_units_area_list(version_units)
    area_list = v_units_area_list['area_list']
    area_list_str = v_units_area_list['area_list_str']



    # load data and init vars, tables, and slices
    md = MetaData()
    db = md.db_base_loader(['units'])
    units = db['units']

    src_filename = md.proc_dest_path(path.join('BehavioralUnits', 'PCA', version_pca,
                                               behunit_params_str(version_fr, timebin, timestep, t_start, t_end),
                                               pca_filter_params_str(filter_params_list_str, counts_thr, area_list_str), 'filter'),
                                     'conditions_events_filter_dict.pkl')
    target_filename = md.proc_dest_path(path.join('BehavioralUnits', 'PCA', version_pca,
                                                  behunit_params_str(version_fr, timebin, timestep, t_start, t_end),
                                                  pca_filter_params_str(filter_params_list_str, counts_thr, area_list_str), 'results'),
                                        'pca_results.pkl')
    print(target_filename)
    if path.exists(target_filename):
        exit()

    conditions_events_filter_dict = md.np_loader(src_filename)

    X = conditions_events_filter_dict['pca_matrix'].transpose()
    X_s = StandardScaler().fit_transform(X.transpose()).transpose()


    # num_conditions = conditions_events_filter_dict['unit_condition_fr_dims'][0]
    # num_timepoints = conditions_events_filter_dict['unit_condition_fr_dims'][1]
    # num_units = len(conditions_events_filter_dict['unit_ind_list'])
    # X_pre_dpca = np.array([np.ravel(unit_array.reshape(num_conditions, num_timepoints).transpose())
    #                        for unit_array
    #                        in X.transpose()]).reshape(num_units, num_timepoints, num_conditions)
    # dpca = dPCA.dPCA(labels='st', regularizer='auto')


    # unit_condition_fr_dims = conditions_events_filter_dict['unit_condition_fr_dims']
    # unit_ind_list = conditions_events_filter_dict['unit_ind_list']
    pca = PCA()
    X_pca = pca.fit_transform(X).transpose()
    X_s_pca = pca.fit_transform(X_s).transpose()

    pca_results = {'X': X, 'X_pca': X_pca,
                   'X_s': X_s, 'X_s_pca': X_s_pca}

    md.np_saver(pca_results, target_filename)



main()
