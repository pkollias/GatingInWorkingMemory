import sys
from rec_analyses import *


def main():
    """ class, balance, fr, counts_thr, area_list, subject, sess_ratio, units_ratio, [overwrite] """
    args_version = sys.argv[1:]
    # args_version = ['job_id=0', 'overwrite=True']
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
    if version['class'] == 'GatingPreBool':
        train_inds_list = [train_inds]
        test_inds_list = [test_inds]
    # split to groups if training on memory or sensory
    elif version['class'] == 'GatedStimulus' or version['class'] == 'Stimulus':
        units_events_valid['Group'] = units_events_valid.apply(lambda row: row[assembly_condition_columns[0]][1], axis=1)
        def assign_index_to_list(series, list):
            for ind, group in series.iteritems():
                list[int(group) - 1].append(ind)
        train_inds_list = [[], []]
        test_inds_list = [[], []]
        assign_index_to_list(units_events_valid.loc[train_inds]['Group'], train_inds_list)
        assign_index_to_list(units_events_valid.loc[test_inds]['Group'], test_inds_list)

    for var, fn in zip([(X_inds_array, y), [train_inds_list, test_inds_list]], target_filename):
        md.np_saver(var, fn)


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
            args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr,
                                                                 args_area_list, args_subject, args_sess_ratio, args_units_ratio)))))

    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))

    return args_version_from_job


main()
