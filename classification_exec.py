import sys
from rec_analyses import *
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis


def main():

    args_version = sys.argv[1:]

    """ class=, balance, fr=, counts_thr=, area_list=, subject=, area=, mode=, mode_seed=, pseudo_seed=, split=, split_ind=, [overwrite=] """
    # args_version = ['class=Stimulus', 'balance=StageGatingCentered', 'fr=ConcatFactor2', 'counts_thr=15',
    #                 'area_list=PFC_Stri', 'subject=Gonzo_Oscar', 'area=PFC', 'mode=Normal',
    #                 'mode_seed=0', 'pseudo_seed=0', 'split=StratifiedStim', 'split_ind=0']
    # args_version = ['job_id=0', 'overwrite=True']

    # load analysis parameters
    version = job_scheduler(args_version, args_from_parse_func)

    # create analysis object
    classifier = ClassificationAnalysis(DataBase([]), version)
    db, md = classifier.db, classifier.db.md

    # overwrite check
    target_filename = classifier.get_path_base('cv_score', classifier.get_train_test_stem())
    print(target_filename)
    if path.exists(target_filename) and ('overwrite' not in version.keys() or not eval(version['overwrite'])):
        exit()

    # init params
    X_inds_array, y = md.np_loader(classifier.get_path_base('X_y', classifier.get_pseudo_session_stem()))
    train_test = md.np_loader(classifier.get_path_base('train_test', classifier.get_train_test_stem()))
    if type(train_test) == dict:
        train_test = train_test[list(train_test.keys())[int(classifier.version['split_ind'])]]
    shrinkage = 'auto'
    solver = 'eigen'


    # load pbt and crop timeseries
    pbt = md.np_loader(classifier.get_path_base('pbt', classifier.get_assemble_stem()))
    if version['fr'] == 'ConcatFactor2':
        pbt.crop_timeseries(-50, 1000)
    t = pbt.timebin_interval.split_to_bins_offset()
    n_t = pbt.timebin_interval.num_of_bins()

    score = np.empty((n_t, n_t, len(train_test)))
    estimator_list = []
    estimator_timeseries = ObjectTimeseries(pbt.timebin_interval)
    # for every timepoint build a training set of classifiers
    for t_train_ind, t_train_i in enumerate(t):

        # for every split
        for split_ind, split_i in enumerate(train_test):

            train_inds = split_i[0]
            X_train = pbt.to_pseudosession_firing_rate(X_inds_array, t_train_ind)[train_inds, :]
            y_train = y[train_inds]
            model = LinearDiscriminantAnalysis(solver=solver, shrinkage=shrinkage)
            model.fit(X_train, y_train)
            estimator_list.append(model)

            # test against for classification score against all other timepoints
            for t_test_ind, t_test_i in enumerate(t):

                test_inds = split_i[1]
                X_test = pbt.to_pseudosession_firing_rate(X_inds_array, t_test_ind)[test_inds, :]
                y_test = y[test_inds]
                score[t_train_ind, t_test_ind, split_ind] = model.score(X_test, y_test)

    md.np_saver(score, classifier.get_path_base('score', classifier.get_train_test_stem()))

    cv_score = np.mean(score, axis=2)
    md.np_saver(cv_score, target_filename)

    estimator_timeseries.set_series(estimator_list)
    md.np_saver(estimator_timeseries, classifier.get_path_base('estimator', classifier.get_train_test_stem()))


def args_from_parse_func(parse_version):

    args_version_list = []



    # for area_list, area in [('PFC', 'PFC'), ('Stri', 'Stri')]:
    #     for session in range(42):
    #         args_class = ['class=GatingPreBool']
    #         args_balance = ['balance=Stimulus']
    #         args_fr = ['fr=WindowGatingClassify', 'fr=ConcatFactor2']
    #         args_counts_thr = ['counts_thr=15']
    #         args_area_list = ['area_list={0:s}'.format(area_list)]
    #         args_subject = ['subject={0:d}'.format(session)]
    #         args_area = ['area={0:s}'.format(area)]
    #         args_mode = ['mode=Normal']
    #         args_mode_seed = ['mode_seed=0']
    #         args_pseudo_seed = ['pseudo_seed=0']
    #         args_split = ['split=StratifiedStim']
    #         args_split_ind = ['split_ind=0']
    #         args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr,
    #                                                              args_area_list, args_subject, args_area, args_mode,
    #                                                              args_mode_seed, args_pseudo_seed, args_split, args_split_ind)))))
    #
    # for balance in ['StageGatingCenteredGatingOnly', 'StageGatingCenteredPostDist1Only']:
    #     for session in range(42):
    #         args_class = ['class=GatedStimulus']
    #         args_balance = ['balance={0:s}'.format(balance)]
    #         args_fr = ['fr=WindowMemoryClassify', 'fr=ConcatFactor2']
    #         args_counts_thr = ['counts_thr=15']
    #         args_area_list = ['area_list=PFC']
    #         args_subject = ['subject={0:d}'.format(session)]
    #         args_area = ['area=PFC']
    #         args_mode = ['mode=Normal']
    #         args_mode_seed = ['mode_seed=0']
    #         args_pseudo_seed = ['pseudo_seed=0']
    #         args_split = ['split=StratifiedBalanceSplit']
    #         args_split_ind = ['split_ind=0']
    #         args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr,
    #                                                              args_area_list, args_subject, args_area, args_mode,
    #                                                              args_mode_seed, args_pseudo_seed, args_split, args_split_ind)))))

    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))

    return args_version_from_job


main()
