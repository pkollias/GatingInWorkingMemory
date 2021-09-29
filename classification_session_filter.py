import sys
from rec_analyses import *


def main():

    args_version = sys.argv[1:]

    """ class=, balance, fr=, counts_thr=, area_list=, subject=, sess_ratio=, units_ratio=, [overwrite=] """
    # args_version = ['class=Stimulus', 'balance=StageGatingCentered', 'fr=ConcatFactor2', 'counts_thr=30',
    # 'area_list=PFC_Stri', 'subject=Gonzo_Oscar', 'sess_ratio=0.5', 'units_ratio=0.5']
    # args_version = ['job_id=0']

    # load analysis parameters
    version = job_scheduler(args_version, args_from_parse_func)

    # create analysis object
    classifier = ClassificationAnalysis(DataBase(['sessions', 'units', 'events', 'conditions']), version)
    db, md = classifier.db, classifier.db.md

    # overwrite check
    target_filename = classifier.get_path_base('filter', classifier.get_filter_session_stem())
    print(target_filename)
    if path.exists(target_filename) and ('overwrite' not in version.keys() or not eval(version['overwrite'])):
        exit()

    # create unit selection filters
    units = db.tables['units']
    # valid units
    valid_units_behavioral_lists = md.np_loader(classifier.get_path_base('valid_units', classifier.get_wrangle_stem()))
    valid_units = valid_units_behavioral_lists.apply(bool)
    # single units
    single_units = units['UnitNum'].ne(0) & units['RatingCode'].ne(7)
    # area list units
    area_list = classifier.version['area_list'].split('_')
    area_list_units = units['Area'].isin(area_list)
    # subject units
    sessions = db.tables['sessions']
    session_entry = sessions.iloc[int(classifier.version['subject'])]
    subject_sessions = [session_entry.name]
    sess = session_entry['Session']
    subject_units = units['Session'].isin(subject_sessions)
    # apply filters
    filter_mask = valid_units & single_units & area_list_units & subject_units
    # get filtered units
    units_sess = units.loc[filter_mask]

    # get set of events that filtered units observe
    events_set = set.union(*list(valid_units_behavioral_lists.loc[units_sess.index].apply(set)))
    # filter by ratio of observation
    sess_ratio = float(version['sess_ratio'])
    ratio_mask = valid_units_behavioral_lists.apply(lambda x: np.mean([el in x for el in events_set]) >= sess_ratio)
    # get list of filtered unit indices
    units_index = units.loc[filter_mask & ratio_mask].index

    # create units_events table
    units_events_filter = db.get_units_events()[sess].loc[sorted(events_set)][sorted(units_index)]
    # create column with valid events
    units_events_filter.columns = units_events_filter.columns.to_list()
    num_units = len(units_events_filter.columns)
    num_units_threshold = float(version['units_ratio']) * num_units
    units_events_filter['valid'] = units_events_filter.sum(axis=1).gt(num_units_threshold)

    md.np_saver(units_events_filter, target_filename)


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
