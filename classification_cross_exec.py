import sys
from rec_analyses import *
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

def main():

    args_version = sys.argv[1:]

    """ class_list=, balance_list=, fr=, counts_thr=, area_list=, subject=, area=, mode=, mode_seed=, pseudo_seed=, split=, 
    classifier_ind=, split_split_ind=, [overwrite=] """
    # args_version = ['class_list=Stimulus_GatedStimulus', 'balance_list=StageGatingPrePost_StageGatingPrePost',
    #                 'fr=ConcatFactor2', 'counts_thr=15', 'area_list=PFC_Stri', 'subject=Gonzo_Oscar', 'area=PFC',
    #                 'mode=Normal', 'mode_seed=0', 'pseudo_seed=0', 'split=StratifiedStim',
    #                 'classifier_ind=0_0', 'split_split_ind=0_0']
    # args_version = ['job_id=0', 'overwrite=True']

    # load analysis parameters
    version = job_scheduler(args_version, args_from_parse_func)

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
    train_key = list(train_dict.keys())[split_split_ind[0]]
    train = train_dict[train_key]

    # get testing classifier
    class_test = classifier_list[classifier_ind[1]]
    db_test, md_test = class_test.db, class_test.db.md
    X_inds_array_train_test_test, y_train_test_test = md_test.np_loader(class_test.get_path_base('X_y', class_test.get_pseudo_session_stem()))
    test_dict = md_test.np_loader(class_test.get_path_base('train_test', class_test.get_train_test_stem()))
    # get split of splits
    test_key = list(test_dict.keys())[split_split_ind[1]]
    test = test_dict[test_key]

    # overwrite check
    stem = (*classifier_list[0].get_train_test_stem(),
            'intermediate',
            '_'.join([class_train.version['class'], train_key, class_test.version['class'], test_key]))
    target_filename = classifier_list[0].get_path_base('cv_score', stem, cross=True)
    print(target_filename, file=sys.stdout)
    if path.exists(target_filename) and ('overwrite' not in version.keys() or not eval(version['overwrite'])):
        exit()

    # load pbt and crop timeseries
    pbt_train = md_train.np_loader(class_train.get_path_base('pbt', class_train.get_assemble_stem()))
    pbt_train.crop_timeseries(-50, 1000)
    pbt_test = md_test.np_loader(class_test.get_path_base('pbt', class_test.get_assemble_stem()))
    pbt_test.crop_timeseries(-50, 1000)

    # get time parameters
    t = pbt_train.timebin_interval.split_to_bins_offset()
    n_t = pbt_train.timebin_interval.num_of_bins()
    # for every timepoint build a training set of classifiers

    score = np.empty((n_t, n_t, len(train)))
    # for every split
    for split_ind, (split_train_i, split_test_i) in enumerate(zip(train, test)):

        # for every timepoint build a set of training classifiers
        for t_train_ind, t_train_i in enumerate(t):

            train_inds = split_train_i[0]
            X_train = pbt_train.to_pseudosession_firing_rate(X_inds_array_train_test_train, t_train_ind)[train_inds, :]
            y_train = y_train_test_train[train_inds]
            model = LinearDiscriminantAnalysis(solver=solver, shrinkage=shrinkage)
            model.fit(X_train, y_train)

            # for every timepoint test against set of classifiers
            for t_test_ind, t_test_i in enumerate(t):

                test_inds = split_test_i[1]
                X_test = pbt_test.to_pseudosession_firing_rate(X_inds_array_train_test_test, t_test_ind)[test_inds, :]
                y_test = y_train_test_test[test_inds]
                score[t_train_ind, t_test_ind, split_ind] = model.score(X_test, y_test)

    cv_score = np.mean(score, axis=2)
    md_train.np_saver(cv_score, target_filename)


def args_from_parse_func(parse_version):

    args_version_list = []

    for args_class_list, args_balance_list in [(['class_list=Stimulus_GatedStimulus'], ['balance_list=StageGatingPrePostMemory_StageGatingPrePostSensory'])]:
        for counts_thr in ['12']:  # ['15', '12', '9']:
            for area_list in ['PFC', 'Stri', 'IT']:  # ['PFC', 'Stri', 'IT']:
                for area in area_list.split('_'):
                    args_fr = ['fr=ConcatFactor2']
                    args_counts_thr = ['counts_thr={0:s}'.format(counts_thr)]
                    args_area_list = ['area_list={0:s}'.format(area_list)]
                    args_subject = ['subject=Gonzo_Oscar']
                    args_area = ['area={0:s}'.format(area)]
                    args_mode = ['mode=Normal']
                    args_mode_seed = ['mode_seed=0']
                    args_pseudo_seed = ['pseudo_seed={0:d}'.format(ps_i) for ps_i in range(5)]
                    args_split = ['split=StratifiedBalanceSplit']  # ['split=StratifiedBalanceSplit', 'split=StratifiedBalanceSplit_StimHalf']
                    args_classifier_ind = ['classifier_ind={0:d}_{1:d}'.format(*t) for t in list(product(range(2), range(2)))]
                    args_split_split_ind = ['split_split_ind={0:d}_{1:d}'.format(*t) for t in list(product(range(3), range(3)))]
                    args_version_list.extend(list(map(list, list(product(args_class_list, args_balance_list, args_fr, args_counts_thr,
                                                                         args_area_list, args_subject, args_area,
                                                                         args_mode, args_mode_seed, args_split,
                                                                         args_classifier_ind, args_split_split_ind, args_pseudo_seed)))))

    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))

    return args_version_from_job


main()
