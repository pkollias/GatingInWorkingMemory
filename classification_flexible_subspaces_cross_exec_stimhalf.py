import sys
from rec_analyses import *
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

def main():
    """ class_list, balance_list, fr, counts_thr, area_list, subject, area, mode, mode_seed, pseudo_seed, split,
    classifier_ind, split_split_ind, shuffle, time_t, [overwrite] """
    args_version = sys.argv[1:]
    # args_version = ['job_id=0', 'overwrite=True']
    version = job_scheduler(args_version, args_from_parse_func)

    classifier_ind = list(map(int, version['classifier_ind'].split('_')))
    split_split_ind = list(map(lambda x: [int(el) for el in list(x)], version['split_split_ind'].split('_')))
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
    train_key_abstract = itemgetter(*split_split_ind[0])(list(np.unique(list(map(itemgetter(0), train_dict.keys())))))  ### 2GROUP
    train_key_abstract = list(train_key_abstract) if type(train_key_abstract) is tuple else [train_key_abstract]
    train_key_combs_list = list(product([train_key_abstract], [[1], [2]]))  ### 2GROUP

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
    test_key_abstract = itemgetter(*split_split_ind[1])(list(np.unique(list(map(itemgetter(0), test_dict.keys())))))  ### 2GROUP
    test_key_abstract = list(test_key_abstract) if type(test_key_abstract) is tuple else [test_key_abstract]
    test_key_combs_list = list(product([test_key_abstract], [[1], [2]]))  ### 2GROUP

    # overwrite check
    stem = (*classifier_list[0].get_train_test_stem(),
            'intermediate',
            '_'.join([class_train.version['class'], ''.join(sorted(train_key_abstract)), class_test.version['class'], ''.join(sorted(test_key_abstract))]))  ### 2GROUP
    fname_stem = '_{0:s}'.format(version['train_t']) if version['train_t'] in ['crossbin', 'sens', 'mem', 'diagonal'] else ''
    fname = ('cv_score' if not int(version['shuffle']) else 'cv_score_{0:04d}'.format(int(version['shuffle']))) + fname_stem
    target_filename = classifier_list[0].get_path_base(fname, stem, cross=True)
    print(target_filename, file=sys.stdout)
    if path.exists(target_filename) and ('overwrite' not in version.keys() or not eval(version['overwrite'])):
        exit()

    # load pbt and crop timeseries
    pbt_train = md_train.np_loader(class_train.get_path_base('pbt', class_train.get_assemble_stem()))
    pbt_test = md_test.np_loader(class_test.get_path_base('pbt', class_test.get_assemble_stem()))

    pbt_train.crop_timeseries(-250, 1550)
    pbt_test.crop_timeseries(-250, 1550)


    # get time parameters
    t = pbt_train.timebin_interval.split_to_bins_offset()


    t_train_list_dict = {'crossbin': [list(range(250, 850, 150))],
                         'sens': [list(range(250, 250 + (2 + 1) * 150, 150))],
                         'mem': [list(range(750, 750 + (2 + 1) * 150, 150))],
                         'diagonal': [[t_i] for t_i in pbt_train.timebin_interval.split_to_bins_offset()],
                         'all': [[t_i] for t_i in pbt_train.timebin_interval.split_to_bins_offset()]}
    t_train_list = t_train_list_dict.get(version['train_t'], [[t_i] for t_i in pbt_train.timebin_interval.split_to_bins_offset()])
    t_inds_train_list = [[t.index(t_i) for t_i in t_list_i] for t_list_i in t_train_list]
    n_t_train = len(t_train_list)

    if version['train_t'] == 'diagonal':
        t_test_list = None
        t_inds_test_list = None
        n_t_test = 1
    else:
        t_test_list = [[t_i] for t_i in pbt_test.timebin_interval.split_to_bins_offset()]
        t_inds_test_list = [[t.index(t_i) for t_i in t_list_i] for t_list_i in t_test_list]
        n_t_test = len(t_test_list)


    # get number of train_test splits
    n_splits = len(list(train_dict.values())[0])

    # for every timepoint build a training set of classifiers
    score = np.empty((n_t_train, n_t_test, n_splits, 2))  ### 2GROUP

    # for all iterations of the keys
    for train_test_ii, (train_key_combs, test_key_combs) in enumerate(zip(train_key_combs_list, test_key_combs_list)):  ### 2GROUP

        train_key_list = list(product(*train_key_combs))
        test_key_list = list(product(*test_key_combs))

        # for every split
        for split_ind in range(n_splits):

            # for every timepoint build a set of training classifiers
            for train_i, (t_train, t_inds_train) in enumerate(zip(t_train_list, t_inds_train_list)):

                X_train_list = []
                y_train_list = []
                for t_train_ind, t_train_i in zip(t_inds_train, t_train):
                    split_train_getter = 0

                    # def balancing_seed():
                    #     return get_seed((version, train_key, split_ind))

                    # balanced sampling
                    train_inds = np.concatenate([train_dict[train_key][split_ind][split_train_getter] for train_key in train_key_list])
                    X_train_list.append(pbt_train.to_pseudosession_firing_rate(X_inds_array_train_test_train, t_train_ind)[train_inds, :])
                    y_train_list.append(y_train_test_train[train_inds])

                X_train = np.vstack(X_train_list)
                y_train = np.hstack(y_train_list)
                model = LinearDiscriminantAnalysis(solver=solver, shrinkage=shrinkage)
                model.fit(X_train, y_train)

                # Check if we are creating the diagonal
                if version['train_t'] == 'diagonal':
                    t_test_list_spec = [t_train]
                    t_inds_test_list_spec = [t_inds_train]
                else:
                    t_test_list_spec = t_test_list
                    t_inds_test_list_spec = t_inds_test_list

                for test_i, (t_test, t_inds_test) in enumerate(zip(t_test_list_spec, t_inds_test_list_spec)):

                    X_test_list = []
                    y_test_list = []
                    for t_test_ind, t_test_i in zip(t_inds_test, t_test):
                        split_test_getter = 1
                        test_inds = np.concatenate([test_dict[test_key][split_ind][split_test_getter] for test_key in test_key_list])
                        X_test_list.append(pbt_test.to_pseudosession_firing_rate(X_inds_array_train_test_test, t_test_ind)[test_inds, :])
                        y_test_list.append(y_train_test_test[test_inds])

                    X_test = np.vstack(X_test_list)
                    y_test = np.hstack(y_test_list)
                    score[train_i, test_i, split_ind, train_test_ii] = model.score(X_test, y_test)  ### 2GROUP

    cv_score = np.mean(score, axis=(2, 3))
    md_train.np_saver(cv_score, target_filename)


def args_from_parse_func(parse_version):

    args_version_list = []

    for args_class_list, args_balance_list in [(['class_list=Stimulus_GatedStimulus'], ['balance_list=StageGatingPrePostMemory_StageGatingPrePostSensory'])]:
        for counts_thr in ['12']:  # ['15', '12', '9']:
            for classifier_ind in [(0, 0)]:
                for train_t, split_split_ind in [('sens', ('012', '0')), ('sens', ('012', '1')), ('sens', ('012', '2')),
                                                 ('mem', ('12', '0')), ('mem', ('12', '1')), ('mem', ('12', '2')),
                                                 ('mem', ('0', '0')), ('mem', ('0', '1')), ('mem', ('0', '2')),
                                                 ('diagonal', ('0', '0')), ('diagonal', ('1', '1')), ('diagonal', ('2', '2'))]:
                    for shuffle in range(1): # 501
                        for mode_seed in range(1000):
                            for area_list in ['PFC']:
                                for area in area_list.split('_'):
                                    args_fr = ['fr=ConcatFactor2']
                                    args_counts_thr = ['counts_thr={0:s}'.format(counts_thr)]
                                    args_area_list = ['area_list={0:s}'.format(area_list)]
                                    args_subject = ['subject=Gonzo_Oscar']
                                    args_area = ['area={0:s}'.format(area)]
                                    args_mode = ['mode=Normal']
                                    args_mode_seed = ['mode_seed={0:d}'.format(mode_seed)]
                                    args_pseudo_seed = ['pseudo_seed=0']
                                    args_split = ['split=StratifiedBalanceSplit_StimHalf']  # ['split=StratifiedBalanceSplit', 'split=StratifiedBalanceSplit_StimHalf']
                                    args_classifier_ind = ['classifier_ind={0:d}_{1:d}'.format(*classifier_ind)]
                                    args_split_split_ind = ['split_split_ind={0:s}_{1:s}'.format(*split_split_ind)]
                                    args_shuffle = ['shuffle={0:02d}'.format(shuffle)]
                                    args_train_t = ['train_t={0:s}'.format(train_t)]
                                    args_version_list.extend(list(map(list, list(product(args_class_list, args_balance_list, args_fr, args_counts_thr,
                                                                                         args_area_list, args_subject, args_area,
                                                                                         args_mode, args_mode_seed, args_split,
                                                                                         args_classifier_ind, args_split_split_ind,
                                                                                         args_pseudo_seed, args_shuffle, args_train_t)))))

    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))

    return args_version_from_job


main()
