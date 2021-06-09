import sys
from rec_analyses import *
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis


def main():

    args_version = sys.argv[1:]

    """ class=, balance, fr=, counts_thr=, area_list=, subject=, area=, mode=, mode_seed=, pseudo_seed=, split=, [overwrite=] """
    # args_version = ['class=GatingPreBool', 'balance=Stimulus', 'fr=ConcatFactor2', 'counts_thr=15',
    #                 'area_list=PFC_Stri', 'subject=Gonzo_Oscar', 'area=PFC', 'mode=Normal',
    #                 'mode_seed=0', 'pseudo_seed=0', 'split=StratifiedStim']
    # args_version = ['job_id=0', 'overwrite=True']

    # load analysis parameters
    version = job_scheduler(args_version, args_from_parse_func)

    # create analysis object
    classifier = ClassificationAnalysis(DataBase([]), version)
    db, md = classifier.db, classifier.db.md

    # load population behavioral timeseries
    pbt = md.np_loader(classifier.get_path_base('pbt', classifier.get_assemble_stem()))
    pbt.crop_timeseries(-50, 1000)
    timebin_interval = pbt.timebin_interval
    t = timebin_interval.split_to_bins_offset()
    pseudosession_inds_array, class_array = classifier.generate_pseudoarray_inds(pbt, int(version['pseudo_seed']))

    X_inds_array = pseudosession_inds_array
    y = class_array
    train_test = classifier.get_train_test_splits_list(pbt)

    md.np_saver((X_inds_array, y), classifier.get_path_base('X_y', classifier.get_pseudo_session_stem()))
    md.np_saver(train_test, classifier.get_path_base('train_test', classifier.get_train_test_stem()))


def args_from_parse_func(parse_version):

    args_class = ['class=GatingPreBool']
    args_balance = ['balance=Stimulus']
    args_fr = ['fr=ConcatFactor2']
    args_counts_thr = ['counts_thr=15']
    args_area_list = ['area_list=PFC_Stri']
    args_subject = ['subject=Gonzo_Oscar']
    args_area = ['area=PFC', 'area=Stri']
    args_mode = ['mode=Bootstrap', 'mode=BootstrapEvents']
    args_mode_seed = ['mode_seed={0:d}'.format(ii) for ii in range(100)]
    args_pseudo_seed = ['pseudo_seed=0']
    args_split = ['split=StratifiedStim', 'split=OneStimTest', 'split=OneStimTrain',
                  'split=WithinGroupTransfer', 'split=AcrossGroupTransfer']
    args_version_list = list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr,
                                                    args_area_list, args_subject, args_area, args_mode,
                                                    args_mode_seed, args_pseudo_seed, args_split))))
    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))

    return args_version_from_job


main()
