import sys
from rec_utils import *
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn import svm
from sklearn.preprocessing import StandardScaler
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer

def main():

    args_version = sys.argv[1:]

    """ class=, balance, fr=, counts_thr=, area_list=, subject=, sess_ratio=, units_ratio=, scaler=, imputer=, classifier=, [overwrite=] """
    # args_version = ['class=Stimulus', 'balance=StageGatingCentered', 'fr=ConcatFactor2', 'counts_thr=30',
    # 'area_list=PFC_Stri', 'subject=Gonzo_Oscar', 'sess_ratio=0.5', 'units_ratio=0.5',
    # 'scaler=standard', 'imputer=multiple', 'classifier=lda']
    # args_version = ['job_id=0', 'overwrite=True']

    # load analysis parameters
    version = job_scheduler(args_version, args_from_parse_func)

    # create analysis object
    classifier = ClassificationAnalysis(DataBase(['units']), version)
    db, md = classifier.db, classifier.db.md

    # overwrite check
    target_filename = classifier.get_path_base('test_result', classifier.get_train_test_session_stem())
    print(target_filename)
    if path.exists(target_filename) and ('overwrite' not in version.keys() or not eval(version['overwrite'])):
        exit()

    # load variables from db
    units_events_filter = md.np_loader(classifier.get_path_base('filter', classifier.get_filter_session_stem()))
    X_inds_array, y = md.np_loader(classifier.get_path_base('X_y', classifier.get_session_stem()))
    [train_inds, test_inds] = md.np_loader(classifier.get_path_base('train_test', classifier.get_train_test_session_stem()))

    # params init
    timebin_interval = timebin_interval_from_version_fr(version['fr'])
    t = timebin_interval.split_to_bins_offset()

    # classifier params
    if version['classifier'] == 'lda':
        classifier_params = {'shrinkage': 'auto',
                             'solver': 'eigen'}
        model_class = LinearDiscriminantAnalysis
    elif version['classifier'] == 'svm':
        classifier_params = {'kernel': 'linear',
                             'probability': True}
        model_class = svm.SVC

    test_result_list = []
    test_result_timeseries = ObjectTimeseries(timebin_interval)
    # for every timepoint build a training set of classifiers
    for t_ind, t_i in enumerate(t):

        # transform (scale, impute) based on params
        units_events_fr = fr_from_units_events(units_events_filter, version['fr'], db, t_ind)
        X_full = units_events_fr.to_numpy()
        # scale
        if version['scaler'] == 'standard':
            X_full = StandardScaler().fit_transform(X_full)
        # impute
        if version['imputer'] == 'multiple':
            imp = IterativeImputer(max_iter=100, random_state=0)
            X_impute = imp.fit_transform(X_full)
        # select valid
        X = X_impute[list(units_events_filter['valid']), :]

        X_train, y_train = X[train_inds, :], y[train_inds]
        X_test, y_test = X[test_inds, :], y[test_inds]
        model = model_class(**classifier_params)
        model.fit(X_train, y_train)

        score = model.score(X_test, y_test)
        predict = model.predict(X_test)
        proba = model.predict_proba(X_test)

        test_result_list.append((score, predict, proba))

    # save estimator
    test_result_timeseries.set_series(test_result_list)
    md.np_saver(test_result_timeseries, target_filename)


def args_from_parse_func(parse_version):

    args_version_list = []

    for session in range(42):
        args_class = ['class=GatingPreBool']
        args_balance = ['balance=Stimulus']
        args_fr = ['fr=WindowGatingClassify']
        args_counts_thr = ['counts_thr=6', 'counts_thr=9', 'counts_thr=15']
        args_area_list = ['area_list=PFC', 'area_list=Stri']
        args_subject = ['subject={0:d}'.format(session)]
        args_sess_ratio = ['sess_ratio={0:s}'.format(ratio) for ratio in ['0', '0.25', '0.5', '0.75', '0.9']]
        args_units_ratio = ['units_ratio={0:s}'.format(ratio) for ratio in ['0', '0.5', '0.6', '0.7', '0.8', '0.9']]
        args_scaler = ['scaler=standard', 'scaler=none']
        args_imputer = ['imputer=multiple']
        args_classifier = ['classifier=lda', 'classifier=svm']
        args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr,
                                                             args_area_list, args_subject, args_sess_ratio, args_units_ratio,
                                                             args_scaler, args_imputer, args_classifier)))))

    for session in range(42):
        args_class = ['class=GatedStimulus']
        args_balance = ['balance=StageGatingCenteredSensoryGatingOnly', 'balance=StageGatingCenteredSensoryPostDist1Only']
        args_fr = ['fr=WindowMemoryClassify', 'fr=WindowInterferenceClassify']
        args_counts_thr = ['counts_thr=6', 'counts_thr=9', 'counts_thr=15']
        args_area_list = ['area_list=PFC']
        args_subject = ['subject={0:d}'.format(session)]
        args_sess_ratio = ['sess_ratio={0:s}'.format(ratio) for ratio in ['0', '0.25', '0.5', '0.75', '0.9']]
        args_units_ratio = ['units_ratio={0:s}'.format(ratio) for ratio in ['0', '0.5', '0.6', '0.7', '0.8', '0.9']]
        args_scaler = ['scaler=standard', 'scaler=none']
        args_imputer = ['imputer=multiple']
        args_classifier = ['classifier=lda', 'classifier=svm']
        args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr,
                                                             args_area_list, args_subject, args_sess_ratio, args_units_ratio,
                                                             args_scaler, args_imputer, args_classifier)))))

    for session in range(42):
        args_class = ['class=Stimulus']
        args_balance = ['balance=StageGatingCenteredMemoryPostDist1Only']
        args_fr = ['fr=WindowMemoryClassify', 'fr=WindowInterferenceClassify']
        args_counts_thr = ['counts_thr=6', 'counts_thr=9', 'counts_thr=15']
        args_area_list = ['area_list=PFC']
        args_subject = ['subject={0:d}'.format(session)]
        args_sess_ratio = ['sess_ratio={0:s}'.format(ratio) for ratio in ['0', '0.25', '0.5', '0.75', '0.9']]
        args_units_ratio = ['units_ratio={0:s}'.format(ratio) for ratio in ['0', '0.5', '0.6', '0.7', '0.8', '0.9']]
        args_scaler = ['scaler=standard', 'scaler=none']
        args_imputer = ['imputer=multiple']
        args_classifier = ['classifier=lda', 'classifier=svm']
        args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr,
                                                             args_area_list, args_subject, args_sess_ratio, args_units_ratio,
                                                             args_scaler, args_imputer, args_classifier)))))

    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))

    return args_version_from_job


main()
