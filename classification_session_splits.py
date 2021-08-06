import sys
from rec_analyses import *


def main():

    args_version = sys.argv[1:]

    """ class=, balance, fr=, counts_thr=, area_list=, subject=, sess_ratio=, units_ratio=, [overwrite=] """
    # args_version = ['class=Stimulus', 'balance=StageGatingCentered', 'fr=ConcatFactor2', 'counts_thr=30',
    # 'area_list=PFC_Stri', 'subject=Gonzo_Oscar', 'sess_ratio=0.5', 'units_ratio=0.5']
    # args_version = ['job_id=0', 'overwrite=True']

    # load analysis parameters
    version = job_scheduler(args_version, args_from_parse_func)

    # create analysis object
    classifier = ClassificationAnalysis(DataBase(['units', 'events', 'conditions']), version)
    db, md = classifier.db, classifier.db.md

    # overwrite check
    target_filename = [classifier.get_path_base('X_y', classifier.get_session_stem()),
                       classifier.get_path_base('train_test', classifier.get_train_test_session_stem())]
    print(target_filename)
    if all([path.exists(fn) for fn in target_filename]) and ('overwrite' not in version.keys() or not eval(version['overwrite'])):
        exit()

    # create unit selection filters
    units_events_filter = md.np_loader(classifier.get_path_base('filter', classifier.get_filter_session_stem()))

    # init params
    events_index = md.preproc_imports['events']['index']
    assembly_condition_columns = classifier.get_assembly_condition_columns()
    # add condition columns
    events_conditions = db.tables['events_conditions'][assembly_condition_columns].reset_index()
    units_events_condition = pd.merge(units_events_filter, events_conditions, on=events_index)
    units_events_valid = units_events_condition.loc[units_events_condition['valid']]

    # create X and y (X is simply range(len(units_events_valid)))
    class_condition_columns = classifier.get_class_condition_columns()
    X_inds_array = units_events_valid.index
    y = zip_columns(units_events_valid, class_condition_columns, 'ClassCondition').astype('category').cat.codes.to_numpy()

    # sample equally assembly condition columns for train inds, rest is test
    train_inds = units_events_valid.groupby(assembly_condition_columns).sample(int(version['counts_thr']), random_state=0).index
    test_inds = units_events_valid.index.difference(train_inds)

    for var, fn in zip([(X_inds_array, y), [train_inds, test_inds]], target_filename):
        md.np_saver(var, fn)


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
        args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr,
                                                             args_area_list, args_subject, args_sess_ratio, args_units_ratio)))))

    for session in range(42):
        args_class = ['class=GatedStimulus']
        args_balance = ['balance=StageGatingCenteredSensoryGatingOnly', 'balance=StageGatingCenteredSensoryPostDist1Only']
        args_fr = ['fr=WindowMemoryClassify', 'fr=WindowInterferenceClassify']
        args_counts_thr = ['counts_thr=6', 'counts_thr=9', 'counts_thr=15']
        args_area_list = ['area_list=PFC']
        args_subject = ['subject={0:d}'.format(session)]
        args_sess_ratio = ['sess_ratio={0:s}'.format(ratio) for ratio in ['0', '0.25', '0.5', '0.75', '0.9']]
        args_units_ratio = ['units_ratio={0:s}'.format(ratio) for ratio in ['0', '0.5', '0.6', '0.7', '0.8', '0.9']]
        args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr,
                                                             args_area_list, args_subject, args_sess_ratio, args_units_ratio)))))

    for session in range(42):
        args_class = ['class=Stimulus']
        args_balance = ['balance=StageGatingCenteredMemoryPostDist1Only']
        args_fr = ['fr=WindowMemoryClassify', 'fr=WindowInterferenceClassify']
        args_counts_thr = ['counts_thr=6', 'counts_thr=9', 'counts_thr=15']
        args_area_list = ['area_list=PFC']
        args_subject = ['subject={0:d}'.format(session)]
        args_sess_ratio = ['sess_ratio={0:s}'.format(ratio) for ratio in ['0', '0.25', '0.5', '0.75', '0.9']]
        args_units_ratio = ['units_ratio={0:s}'.format(ratio) for ratio in ['0', '0.5', '0.6', '0.7', '0.8', '0.9']]
        args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr,
                                                             args_area_list, args_subject, args_sess_ratio, args_units_ratio)))))

    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))

    return args_version_from_job


main()
