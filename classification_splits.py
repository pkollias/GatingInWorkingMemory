import sys
from rec_analyses import *


def main():

    args_version = sys.argv[1:]

    """ class=, balance, fr=, counts_thr=, area_list=, subject=, area=, mode=, mode_seed=, pseudo_seed=, split=, [overwrite=] """
    # args_version = ['class=Stimulus', 'balance=StageGatingCentered', 'fr=ConcatFactor2', 'counts_thr=15',
    #                 'area_list=PFC_Stri', 'subject=Gonzo_Oscar', 'area=PFC', 'mode=Bootstrap',
    #                 'mode_seed=0', 'pseudo_seed=0', 'split=StratifiedStim']
    # args_version = ['job_id=0', 'overwrite=True']

    # load analysis parameters
    version = job_scheduler(args_version, args_from_parse_func)

    # create analysis object
    classifier = ClassificationAnalysis(DataBase([]), version)
    db, md = classifier.db, classifier.db.md

    # load population behavioral timeseries
    pbt = md.np_loader(classifier.get_path_base('pbt', classifier.get_assemble_stem()))
    pseudosession_inds_array, class_array = classifier.generate_pseudoarray_inds(pbt, int(version['pseudo_seed']))

    X_inds_array = pseudosession_inds_array
    y = class_array
    train_test = classifier.get_train_test_splits_list(pbt)

    md.np_saver((X_inds_array, y), classifier.get_path_base('X_y', classifier.get_pseudo_session_stem()))
    md.np_saver(train_test, classifier.get_path_base('train_test', classifier.get_train_test_stem()))


def args_from_parse_func(parse_version):

    args_version_list = []

    for area_list, area in [('PFC', 'PFC'), ('Stri', 'Stri'), ('IT', 'IT')]:
        args_class = ['class=GatingPreBoolGeneralized']
        args_balance = ['balance=Stimulus']
        args_fr = ['fr=ConcatFactor2']
        args_counts_thr = ['counts_thr={0:s}'.format(counts_thr) for counts_thr in ['12', '15']]
        args_area_list = ['area_list={0:s}'.format(area_list)]
        args_subject = ['subject=Gonzo_Oscar']
        args_area = ['area={0:s}'.format(area)]
        args_mode = ['mode=Normal']
        args_mode_seed = ['mode_seed=0']
        # args_mode_seed = ['mode_seed={0:s}'.format(mode_seed) for mode_seed in [str(ms) for ms in range(10)]]
        args_pseudo_seed = ['pseudo_seed={0:s}'.format(pseudo_seed) for pseudo_seed in [str(ps) for ps in range(10)]]
        args_split = ['split={0:s}'.format(split) for split in ['StratifiedStim', 'OneStimTest', 'OneStimTrain', 'WithinGroupTransfer', 'AcrossGroupTransfer']]
        args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr,
                                                             args_area_list, args_subject, args_area, args_mode,
                                                             args_mode_seed, args_pseudo_seed, args_split)))))

    # args_version_list = []
    #
    # for area_list, area in [('PFC', 'PFC'), ('Stri', 'Stri'), ('IT', 'IT')]:
    #     for class_i, balance in [('Stimulus', 'StageGatingPrePostMemory'), ('Stimulus', 'StageGatingCenteredMemory'),
    #                              ('GatedStimulus', 'StageGatingPrePostSensory'), ('GatedStimulus', 'StageGatingCenteredSensory')]:
    #         args_class = ['class={0:s}'.format(class_i)]
    #         args_balance = ['balance={0:s}'.format(balance)]
    #         args_fr = ['fr=ConcatFactor2']
    #         args_counts_thr = ['counts_thr=12']
    #         args_area_list = ['area_list={0:s}'.format(area_list)]
    #         args_subject = ['subject=Gonzo_Oscar']
    #         args_area = ['area={0:s}'.format(area)]
    #         args_mode = ['mode=Normal']
    #         args_mode_seed = ['mode_seed=0']
    #         args_pseudo_seed = ['pseudo_seed={0:d}'.format(ps_i) for ps_i in range(5)]
    #         args_split = ['split=StratifiedBalanceSplit_StimHalf']
    #         args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr,
    #                                                              args_area_list, args_subject, args_area, args_mode,
    #                                                              args_mode_seed, args_pseudo_seed, args_split)))))

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
    #         args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr,
    #                                                              args_area_list, args_subject, args_area, args_mode,
    #                                                              args_mode_seed, args_pseudo_seed, args_split)))))
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
    #         args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr,
    #                                                              args_area_list, args_subject, args_area, args_mode,
    #                                                              args_mode_seed, args_pseudo_seed, args_split)))))

    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))

    return args_version_from_job


main()
