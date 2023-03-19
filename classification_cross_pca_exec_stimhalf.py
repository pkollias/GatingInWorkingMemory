import sys
from rec_analyses import *
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.pipeline import Pipeline
from rec_utils import fbt_df_from_PCA

def main():
    """ class_list, balance_list, fr, counts_thr, area_list, subject, area, mode, mode_seed, pseudo_seed, split,
    classifier_ind, split_split_ind, shuffle, pca_mode, [overwrite] """
    args_version = sys.argv[1:]
    # args_version = ['job_id=0', 'overwrite=True']

    classifier_ind = list(map(int, version['classifier_ind'].split('_')))
    split_split_ind = list(map(int, version['split_split_ind'].split('_')))
    class_list = version['class_list'].split('_')
    balance_list = version['balance_list'].split('_')
    def version_modifier(class_i, balance_i):
        version_mod = version.copy()
        version_mod['class'] = class_i
        version_mod['balance'] = balance_i
        return version_mod
    # create analysis object
    classifier_list = [ClassificationAnalysis(DataBase([]), version_modifier(class_i, balance_i))
                       for class_i, balance_i
                       in zip(class_list, balance_list)]

    shrinkage = 'auto'
    solver = 'eigen'

    # get training classifier
    class_train = classifier_list[classifier_ind[0]]
    db_train, md_train = class_train.db, class_train.db.md
    X_inds_array_train_test_train, y_train_test_train = md_train.np_loader(class_train.get_path_base('X_y', class_train.get_pseudo_session_stem()))
    train_dict = md_train.np_loader(class_train.get_path_base('train_test', class_train.get_train_test_stem()))
    # get split of splits
    train_key_abstract = list(np.unique(list(map(itemgetter(0), train_dict.keys()))))[split_split_ind[0]]  ### 2GROUP
    train_key_list = list(product([train_key_abstract], [1, 2]))  ### 2GROUP

    # if random permutation of labels then shuffle with seed
    if int(version['shuffle']):
        np.random.seed(int(version['shuffle']))
        y_temp = y_train_test_train.reshape((2, -1)).transpose()
        np.random.shuffle(y_temp)
        y_train_test_train = y_temp.transpose().reshape(-1)  ### 2GROUP

    # get testing classifier
    class_test = classifier_list[classifier_ind[1]]
    db_test, md_test = class_test.db, class_test.db.md
    X_inds_array_train_test_test, y_train_test_test = md_test.np_loader(class_test.get_path_base('X_y', class_test.get_pseudo_session_stem()))
    test_dict = md_test.np_loader(class_test.get_path_base('train_test', class_test.get_train_test_stem()))
    # get split of splits
    test_key_abstract = list(np.unique(list(map(itemgetter(0), test_dict.keys()))))[split_split_ind[1]]  ### 2GROUP
    test_key_list = list(product([test_key_abstract], [1, 2]))  ### 2GROUP

    pbt_train_src = md_train.np_loader(class_train.get_path_base('pbt', class_train.get_assemble_stem()))
    pbt_train_src.crop_timeseries(-50, 1000)
    pbt_test_src = md_test.np_loader(class_test.get_path_base('pbt', class_test.get_assemble_stem()))
    pbt_test_src.crop_timeseries(-50, 1000)

    # overwrite check
    # stem = (*classifier_list[0].get_train_test_stem(),
    #         'intermediate',
    #         '_'.join([class_train.version['class'], train_key_abstract, class_test.version['class'], test_key_abstract]),
    #         'pca')  ### 2GROUP
    # fname = 'cv_score' if not int(version['shuffle']) else 'cv_score_{0:04d}'.format(int(version['shuffle']))
    # target_filename = classifier_list[0].get_path_base(fname, stem, cross=True)
    # print(target_filename, file=sys.stdout)
    # if path.exists(target_filename) and ('overwrite' not in version.keys() or not eval(version['overwrite'])):
    #     exit()

    # load pbt and crop timeseries
    # _, pca = md_train.np_loader(class_train.get_path_base('pca', class_train.get_assemble_stem(), cross=True))
    # n_components = next(ind for ind, x in enumerate(pca.explained_variance_ratio_.cumsum()) if x > 2/3)

    # train
    # fname = 'fbt'
    # stem_train = tuple(list(class_train.get_assemble_stem()) + \
    #              ['_'.join([class_train.version[key] for key in ['class', 'balance']])])
    # pbt_train = md_train.np_loader(class_train.get_path_base(fname, stem_train, cross=True))
    # # test
    # stem_test = tuple(list(class_test.get_assemble_stem()) + \
    #              ['_'.join([class_test.version[key] for key in ['class', 'balance']])])
    # pbt_test = md_test.np_loader(class_test.get_path_base(fname, stem_test, cross=True))
    # #reduce dims
    # X_inds_array_train_test_train = X_inds_array_train_test_train[:, :n_components]
    # pbt_train.df = pbt_train.df.loc[pbt_train.df['Factor'].lt(n_components)]
    # X_inds_array_train_test_test = X_inds_array_train_test_test[:, :n_components]
    # pbt_test.df = pbt_test.df.loc[pbt_test.df['Factor'].lt(n_components)]


    # get time parameters
    t = pbt_train_src.timebin_interval.split_to_bins_offset()
    n_t = pbt_train_src.timebin_interval.num_of_bins()
    # for every timepoint build a training set of classifiers

    score = np.empty((n_t, n_t, len(list(train_dict.values())[0]), 2))  ### 2GROUP

    # for all iterations of the keys
    for train_test_ii, (train_key, test_key) in enumerate(zip(train_key_list, test_key_list)):  ### 2GROUP
        train = train_dict[train_key]  ### 2GROUP
        test = test_dict[test_key]  ### 2GROUP

        # for every split
        for split_ind, (split_train_i, split_test_i) in enumerate(zip(train, test)):

            # estimate pca fitting df and X
            df_train_fit = pbt_train_src.df.iloc[sorted(X_inds_array_train_test_train[split_train_i[0], :].reshape(-1))]
            pbt_train_fit = pbt_train_src.init_with_df(df_train_fit)
            if version['pca_mode'] == 'mean':
                X_fit, _ = pbt_train_fit.average_instances(['Unit', 'Unit_Code', 'Condition']).to_PCA_array()
            else:
                X_fit, _ = pbt_train_fit.to_PCA_array()

            # train pca model
            scaler = StandardScaler()
            pca = PCA(n_components=2/3)
            pca.fit(scaler.fit_transform(X_fit))
            pipeline = Pipeline([('scaler', scaler), ('pca', pca)])

            # pipeline transform (scale, pca) training data
            # get df rows of indices that correspond to training data
            pbt_train = pbt_train_src
            # transform
            X_train, records_train = pbt_train.to_PCA_array()
            X_factor_train = pipeline.transform(X_train)
            # create fbt
            fbt_train = FactorBehavioralTimeseries(fbt_df_from_PCA(X_factor_train, records_train,
                                                                   X_factor_train.shape[1], pbt_train.timebin_interval),
                                                   pbt_train.condition_labels, pbt_train.timebin_interval)
            pbt_test = pbt_test_src
            # transform
            X_test, records_test = pbt_test.to_PCA_array()
            X_factor_test = pipeline.transform(X_test)
            # create fbt
            fbt_test = FactorBehavioralTimeseries(fbt_df_from_PCA(X_factor_test, records_test,
                                                                   X_factor_test.shape[1], pbt_test.timebin_interval),
                                                   pbt_test.condition_labels, pbt_test.timebin_interval)

            # housekeeping for number of components and fbt index
            # n_components
            n_components = pca.n_components_
            X_inds_train = X_inds_array_train_test_train[:, :n_components]
            X_inds_test = X_inds_array_train_test_test[:, :n_components]
            # index TODO ???

            # for every timepoint build a set of training classifiers
            for t_train_ind, t_train_i in enumerate(t):

                train_inds = split_train_i[0]
                X_train = fbt_train.to_pseudosession_firing_rate(X_inds_train, t_train_ind)[train_inds, :]
                y_train = y_train_test_train[train_inds]
                model = LinearDiscriminantAnalysis(solver=solver, shrinkage=shrinkage)
                model.fit(X_train, y_train)

                # for every timepoint test against set of classifiers
                for t_test_ind, t_test_i in enumerate(t):

                    print(t_test_ind, t_train_ind, split_ind, train_test_ii)
                    test_inds = split_test_i[1]
                    X_test = fbt_test.to_pseudosession_firing_rate(X_inds_test, t_test_ind)[test_inds, :]
                    y_test = y_train_test_test[test_inds]
                    score[t_train_ind, t_test_ind, split_ind, train_test_ii] = model.score(X_test, y_test)  ### 2GROUP

    cv_score = np.mean(score, axis=(2, 3))
    md_train.np_saver(cv_score, target_filename)


def args_from_parse_func(parse_version):

    args_version_list = []

    for args_class_list, args_balance_list in [(['class_list=Stimulus_GatedStimulus'], ['balance_list=StageGatingPrePostMemory_StageGatingPrePostSensory'])]:
        for counts_thr in ['12']:  # ['15', '12', '9']:
            for clasiifier_ind in product(range(2), range(2)):
                for split_split_ind in product(range(3), range(3)):
                    for shuffle in range(1):
                        for area_list in ['PFC', 'Stri', 'IT']:
                            for area in area_list.split('_'):
                                args_fr = ['fr=ConcatFactor2']
                                args_counts_thr = ['counts_thr={0:s}'.format(counts_thr)]
                                args_area_list = ['area_list={0:s}'.format(area_list)]
                                args_subject = ['subject=Gonzo_Oscar']
                                args_area = ['area={0:s}'.format(area)]
                                args_mode = ['mode=Normal']
                                args_mode_seed = ['mode_seed=0']
                                args_pseudo_seed = ['pseudo_seed=0']
                                args_split = ['split=StratifiedBalanceSplit_StimHalf']  # ['split=StratifiedBalanceSplit', 'split=StratifiedBalanceSplit_StimHalf']
                                args_classifier_ind = ['classifier_ind={0:d}_{1:d}'.format(*clasiifier_ind)]
                                args_split_split_ind = ['split_split_ind={0:d}_{1:d}'.format(*split_split_ind)]
                                args_shuffle = ['shuffle={0:02d}'.format(shuffle)]
                                args_pca_mode = ['pca_mode=mean']
                                args_version_list.extend(list(map(list, list(product(args_class_list, args_balance_list, args_fr, args_counts_thr,
                                                                                     args_area_list, args_subject, args_area,
                                                                                     args_mode, args_mode_seed, args_split,
                                                                                     args_classifier_ind, args_split_split_ind,
                                                                                     args_pseudo_seed, args_shuffle, args_pca_mode)))))

    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))

    return args_version_from_job


main()
