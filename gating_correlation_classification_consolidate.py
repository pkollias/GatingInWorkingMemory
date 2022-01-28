import sys
from rec_utils import *


def main():
    """ class, balance, fr, area_list, subject, [overwrite] """
    args_version = sys.argv[1:]
    # args_version = ['job_id=0', 'overwrite=True']
    version = job_scheduler(args_version, args_from_parse_func)

    version_filename = version.copy()

    result_dict = {}

    for args_suffix in product(['6', '9', '15'], ['0', '0.25', '0.5', '0.75', '0.9'],
                               ['0', '0.5', '0.6', '0.7', '0.8', '0.9'], ['standardn', 'none'],
                               ['multiple'], ['lda', 'svm']):

        counts_thr, sess_ratio, units_ratio, scaler, imputer, classifier = args_suffix
        version['counts_thr'] = counts_thr
        version['sess_ratio'] = sess_ratio
        version['units_ratio'] = units_ratio
        version['scaler'] = scaler
        version['imputer'] = imputer
        version['classifier'] = classifier

        # create analysis object
        classifier = ClassificationAnalysis(DataBase(['units']), version)
        db, md = classifier.db, classifier.db.md

        try:
            events = md.np_loader(classifier.get_path_base('events', classifier.get_train_test_session_stem()))
            result_dict[args_suffix] = (events.mean(), events)
        except OSError as e:
            print('File not found')
            pass

    if 'overwrite' in version.keys():
        version_str = '_'.join(list(version_filename.values())[:-1])
    else:
        version_str = '_'.join(version_filename.values())
    target_filename = MetaData().proc_dest_path(path.join('BehavioralUnits', 'GatingCorrelation',
                                                          'consolidate', version_str),
                                                'result_dict.pkl')
    print(target_filename)
    MetaData().np_saver(result_dict, target_filename)


def args_from_parse_func(parse_version):

    args_version_list = []

    for session in range(42):
        args_class = ['class=GatingPreBool']
        args_balance = ['balance=Stimulus']
        args_fr = ['fr=WindowGatingClassify']
        args_area_list = ['area_list=PFC', 'area_list=Stri']
        args_subject = ['subject={0:d}'.format(session)]
        args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_area_list, args_subject)))))

    for session in range(42):
        args_class = ['class=GatedStimulus']
        args_balance = ['balance=StageGatingCenteredSensoryGatingOnly', 'balance=StageGatingCenteredSensoryPostDist1Only']
        args_fr = ['fr=WindowMemoryClassify', 'fr=WindowInterferenceClassify']
        args_area_list = ['area_list=PFC']
        args_subject = ['subject={0:d}'.format(session)]
        args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_area_list, args_subject)))))

    for session in range(42):
        args_class = ['class=Stimulus']
        args_balance = ['balance=StageGatingCenteredMemoryPostDist1Only']
        args_fr = ['fr=WindowMemoryClassify', 'fr=WindowInterferenceClassify']
        args_area_list = ['area_list=PFC']
        args_subject = ['subject={0:d}'.format(session)]
        args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_area_list, args_subject)))))

    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))

    return args_version_from_job


main()
