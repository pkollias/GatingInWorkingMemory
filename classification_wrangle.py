import sys
from rec_analyses import *


def main():
    """ class, balance, fr, counts_thr, [overwrite] """
    args_version = sys.argv[1:]
    # args_version = ['job_id=0']
    version = job_scheduler(args_version, args_from_parse_func)

    # create analysis object
    classifier = ClassificationAnalysis(DataBase(['trials', 'units', 'events', 'conditions']), version)
    db, md = classifier.db, classifier.db.md

    # overwrite check
    target_filename = classifier.get_path_base('valid_units', classifier.get_wrangle_stem())
    print(target_filename)
    if path.exists(target_filename) and ('overwrite' not in version.keys() or not eval(version['overwrite'])):
        exit()

    # mark units as valid based on number of events and
    db.merge_events_conditions_trials()
    units = db.tables['units']
    valid_units_behavioral_lists = units.apply(lambda row: classifier.assess_unit_events(row.name), axis=1)

    md.np_saver(valid_units_behavioral_lists, target_filename)


def args_from_parse_func(parse_version):

    args_version_list = []

    args_class = ['class=Stimulus']
    args_balance = ['balance=StageGatingPrePostMemory']
    args_fr = ['fr=ConcatFactor2']
    args_counts_thr = ['counts_thr=12']
    args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr)))))

    args_class = ['class=GatedStimulus']
    args_balance = ['balance=StageGatingPrePostSensory']
    args_fr = ['fr=ConcatFactor2']
    args_counts_thr = ['counts_thr=12']
    args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr)))))

    # args_class = ['class=GatingPreBoolGeneralizedCue']
    # args_balance = ['balance=Cue']
    # args_fr = ['fr=ConcatFactor2']
    # args_counts_thr = ['counts_thr=1']
    # args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr)))))
    #
    # args_class = ['class=GatingPreBool']
    # args_balance = ['balance=Stimulus']
    # args_fr = ['fr={0:s}'.format(version_fr) for version_fr in ['200_400', '400_600', '600_800', '800_1000']]
    # args_counts_thr = ['counts_thr=12', 'counts_thr=15']
    # args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr)))))
    #
    # args_class = ['class=GatedStimulus']
    # args_balance = ['balance=StageGatingCenteredSensoryGatingOnly', 'balance=StageGatingCenteredSensoryPostDist1Only']
    # args_fr = ['fr={0:s}'.format(version_fr) for version_fr in ['200_400', '400_600', '600_800', '800_1000']]
    # args_counts_thr = ['counts_thr=12', 'counts_thr=15']
    # args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr)))))
    #
    # args_class = ['class=Stimulus']
    # args_balance = ['balance=StageGatingCenteredMemoryPostDist1Only']
    # args_fr = ['fr={0:s}'.format(version_fr) for version_fr in ['200_400', '400_600', '600_800', '800_1000']]
    # args_counts_thr = ['counts_thr=12', 'counts_thr=15']
    # args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr)))))

    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))

    return args_version_from_job


main()
