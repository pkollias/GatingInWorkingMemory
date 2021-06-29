import sys
from rec_analyses import *


def main():

    args_version = sys.argv[1:]

    """ class=, balance, fr=, counts_thr=, area_list=, subject=, [overwrite=] """
    # args_version = ['class=Stimulus', 'balance=StageGatingCentered', 'fr=ConcatFactor2', 'counts_thr=15',
    # 'area_list=PFC_Stri', 'subject=0']
    # args_version = ['job_id=0']

    # load analysis parameters
    version = job_scheduler(args_version, args_from_parse_func)

    # create analysis object
    classifier = ClassificationAnalysis(DataBase(['sessions', 'units', 'events', 'conditions']), version)
    db, md = classifier.db, classifier.db.md

    # overwrite check
    target_filename = classifier.get_path_base('filter', classifier.get_filter_stem())
    print(target_filename)
    if path.exists(target_filename) and ('overwrite' not in version.keys() or not eval(version['overwrite'])):
        exit()

    # create unit selection filters
    units = db.tables['units']
    # valid units
    valid_units_events = md.np_loader(classifier.get_path_base('valid_units', classifier.get_wrangle_stem()))
    valid_units = valid_units_events.apply(bool)
    # single units
    single_units = units['UnitNum'].ne(0) & units['RatingCode'].ne(7)
    # area list units
    area_list = classifier.version['area_list'].split('_')
    area_list_units = units['Area'].isin(area_list)
    # subject units
    sessions = db.tables['sessions']
    subject_sessions = [sessions.iloc[int(classifier.version['subject'])].name]  ### TODO: only difference, can edit base function
    subject_units = units['Session'].isin(subject_sessions)
    # apply filters
    db.tables['units'] = units.loc[valid_units & single_units & area_list_units & subject_units]

    # create behavioral_units table
    units_events = valid_units_events.loc[db.tables['units'].index]
    behavioral_units = pd.DataFrame([(sess, channum, unitnum, trialnum, stageindex)
                                     for (sess, channum, unitnum), events_list in units_events.iteritems()
                                     for _, trialnum, stageindex in events_list],
                                    columns=md.proc_imports['behavioral_units']['index'])

    md.np_saver(behavioral_units, target_filename)


def args_from_parse_func(parse_version):

    args_version_list = []

    # for session in range(42):
    #     args_class = ['class=GatingPreBool']
    #     args_balance = ['balance=Stimulus']
    #     args_fr = ['fr=WindowGatingClassify']
    #     args_counts_thr = ['counts_thr=15']
    #     args_area_list = ['area_list=PFC', 'area_list=Stri']
    #     args_subject = ['subject={0:d}'.format(session)]
    #     args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr,
    #                                                          args_area_list, args_subject)))))

    for session in range(42):
        args_class = ['class=GatedStimulus']
        args_balance = args_balance = ['balance=StageGatingCenteredMemoryGatingOnly', 'balance=StageGatingCenteredMemoryPostDist1Only']
        # ['balance=StageGatingCenteredGatingOnly', 'balance=StageGatingCenteredPostDist1Only']
        args_fr = ['fr=WindowMemoryClassify']
        args_counts_thr = ['counts_thr=6', 'counts_thr=9', 'counts_thr=15']
        args_area_list = ['area_list=PFC']
        args_subject = ['subject={0:d}'.format(session)]
        args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr,
                                                             args_area_list, args_subject)))))

    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))

    return args_version_from_job


main()
