import sys
from rec_analyses import *


def main():

    args_version = sys.argv[1:]

    """ class=, balance, fr=, counts_thr=, area_list=, subject=, [overwrite=] """
    # args_version = ['class=Stimulus', 'balance=StageGatingCentered', 'fr=ConcatFactor2', 'counts_thr=15',
    # 'area_list=PFC_Stri', 'subject=Gonzo_Oscar']
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
    valid_units_behavioral_lists = md.np_loader(classifier.get_path_base('valid_units', classifier.get_wrangle_stem()))
    valid_units = valid_units_behavioral_lists.apply(bool)
    # single units
    single_units = units['UnitNum'].ne(0) & units['RatingCode'].ne(7)
    # area list units
    area_list = classifier.version['area_list'].split('_')
    area_list_units = units['Area'].isin(area_list)
    # subject units
    sessions = db.tables['sessions']
    subject_sessions = sessions.loc[sessions['Subject'].isin(classifier.version['subject'].split('_'))].index
    subject_units = units['Session'].isin(subject_sessions)
    # apply filters
    db.tables['units'] = units.loc[valid_units & single_units & area_list_units & subject_units]

    # create behavioral_units table
    units_behavioral_lists = valid_units_behavioral_lists.loc[db.tables['units'].index]
    behavioral_units = pd.DataFrame([(sess, channum, unitnum, trialnum, stageindex)
                                     for (sess, channum, unitnum), events_list in units_behavioral_lists.iteritems()
                                     for _, trialnum, stageindex in events_list],
                                    columns=md.proc_imports['behavioral_units']['index'])

    md.np_saver(behavioral_units, target_filename)


def args_from_parse_func(parse_version):

    args_version_list = []

    args_class = ['class=GatingPreBoolGeneralizedCue']
    args_balance = ['balance=Cue']
    args_fr = ['fr=ConcatFactor2']
    args_counts_thr = ['counts_thr={0:s}'.format(counts_thr) for counts_thr in ['1']]
    args_area_list = ['area_list={0:s}'.format(area_list) for area_list in ['PFC', 'Stri', 'IT']]
    args_subject = ['subject=Gonzo_Oscar']
    args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr,
                                                         args_area_list, args_subject)))))

    # args_version_list = []
    #
    # for class_i, balance in [('Stimulus', 'StageGatingPrePostMemory'), ('Stimulus', 'StageGatingCenteredMemory'),
    #                          ('GatedStimulus', 'StageGatingPrePostSensory'), ('GatedStimulus', 'StageGatingCenteredSensory')]:
    #     for counts_thr in ['12']:
    #         for area_list in ['PFC', 'Stri', 'IT']:
    #             args_class = ['class={0:s}'.format(class_i)]
    #             args_balance = ['balance={0:s}'.format(balance)]
    #             args_fr = ['fr=ConcatFactor2']
    #             args_counts_thr = ['counts_thr={0:s}'.format(counts_thr)]
    #             args_area_list = ['area_list={0:s}'.format(area_list)]
    #             args_subject = ['subject=Gonzo_Oscar']
    #             args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr,
    #                                                                  args_area_list, args_subject)))))

    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))

    return args_version_from_job


main()
