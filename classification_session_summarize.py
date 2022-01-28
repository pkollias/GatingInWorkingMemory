import sys
from rec_utils import *


def main():
    """ class, balance, fr, counts_thr, area_list, subject, sess_ratio, units_ratio, scaler, imputer, classifier, [overwrite] """
    args_version = sys.argv[1:]
    # args_version = ['job_id=0', 'overwrite=True']
    version = job_scheduler(args_version, args_from_parse_func)

    # create analysis object
    classifier = ClassificationAnalysis(DataBase(['units']), version)
    db, md = classifier.db, classifier.db.md

    target_filename = classifier.get_path_base('events', classifier.get_train_test_session_stem())
    print(target_filename)
    if path.exists(target_filename) and ('overwrite' not in version.keys() or not eval(version['overwrite'])):
        exit()

    # load necessary variables
    test_result_timeseries = md.np_loader(classifier.get_path_base('test_result', classifier.get_train_test_session_stem()))
    units_events_filter = md.np_loader(classifier.get_path_base('filter', classifier.get_filter_session_stem()))
    _, y = md.np_loader(classifier.get_path_base('X_y', classifier.get_session_stem()))
    [_, test_inds_list] = md.np_loader(classifier.get_path_base('train_test', classifier.get_train_test_session_stem()))
    test_inds = list(chain(*test_inds_list))
    # event indices of tested events
    test_event_inds = units_events_filter.loc[units_events_filter['valid']].index[test_inds]

    # save probability of correct class to series
    correct_log_proba = [np.array(test_result_timeseries.series[0][2])[row, col] for row, col in enumerate(y[test_inds] % 2)]
    test_event_multiindex = pd.MultiIndex.from_tuples(test_event_inds)
    events_series = pd.Series(correct_log_proba, name='CorrectLogProba', index=test_event_multiindex)

    md.np_saver(events_series, target_filename)


def args_from_parse_func(parse_version):

    args_version_list = []

    for session in range(42):
        for class_arg, balance, area_list in [('GatingPreBool', 'Stimulus', 'PFC'),
                                              ('GatingPreBool', 'Stimulus', 'Stri'),
                                              ('GatedStimulus', 'StageGatingCenteredSensoryGatingOnly', 'PFC'),
                                              ('GatedStimulus', 'StageGatingCenteredSensoryPostDist1Only', 'PFC'),
                                              ('Stimulus', 'StageGatingCenteredMemoryPostDist1Only', 'PFC')]:
            args_class = ['class={0:s}'.format(class_arg)]
            args_balance = ['balance={0:s}'.format(balance)]
            args_fr = ['fr={0:s}'.format(version_fr) for version_fr in ['200_400', '400_600', '600_800', '800_1000']]
            args_counts_thr = ['counts_thr=12', 'counts_thr=15']
            args_area_list = ['area_list={0:s}'.format(area_list)]
            args_subject = ['subject={0:d}'.format(session)]
            args_sess_ratio = ['sess_ratio={0:s}'.format(ratio) for ratio in ['0.75']]
            args_units_ratio = ['units_ratio={0:s}'.format(ratio) for ratio in ['0.6']]
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
