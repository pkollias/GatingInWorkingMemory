import sys
from rec_analyses import *
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from scipy.stats import ttest_1samp


def main():

    session = 36

    args_version_dict = {'gating': ['class=GatingPreBool', 'balance=Stimulus', 'fr=WindowGatingClassify', 'counts_thr=15',
                               'area_list=PFC', 'subject={0:d}'.format(session), 'area=PFC', 'mode=Normal',
                               'mode_seed=0', 'pseudo_seed=0', 'split=StratifiedStim', 'split_ind=0'],
                         'memory': ['class=GatedStimulus', 'balance=StageGatingCenteredGatingOnly', 'fr=WindowMemoryClassify', 'counts_thr=15',
                                    'area_list=PFC', 'subject={0:d}'.format(session), 'area=PFC', 'mode=Normal',
                                    'mode_seed=0', 'pseudo_seed=0', 'split=StratifiedBalanceSplit', 'split_ind=0']}

    # load analysis parameters
    version_dict = {key: job_scheduler(args_version) for key, args_version in args_version_dict.items()}

    # create analysis object
    classifier_dict = {key: ClassificationAnalysis(DataBase([]), version) for key, version in version_dict.items()}
    db_dict = {key: classifier.db for key, classifier in classifier_dict.items()}
    md_dict = {key: db.md for key, db in db_dict.items()}

    # overwrite check
    cv_score_dict = {key: md_dict[key].np_loader(classifier.get_path_base('cv_score', classifier.get_train_test_stem())) for key, classifier in classifier_dict.items()}
    score_dict = {key: md_dict[key].np_loader(classifier.get_path_base('score', classifier.get_train_test_stem())) for key, classifier in classifier_dict.items()}
    est_dict = {key: md_dict[key].np_loader(classifier.get_path_base('estimator', classifier.get_train_test_stem())) for key, classifier in classifier_dict.items()}
    behavioral_units_dict = {key: md_dict[key].np_loader(classifier.get_path_base('filter', classifier.get_filter_stem())) for key, classifier in classifier_dict.items()}
    pbt_dict = {key: md_dict[key].np_loader(classifier.get_path_base('pbt', classifier.get_assemble_stem())) for key, classifier in classifier_dict.items()}
    units_dict = {key: pbt.get_unit_inds() for key, pbt in pbt_dict.items()}

    db = DataBase(['sessions', 'units', 'events', 'conditions'])
    md = db.md


def args_from_parse_func(parse_version):

    args_version_list = []

    for area_list, area in [('PFC', 'PFC'), ('Stri', 'Stri')]:
        for session in range(42):
            args_class = ['class=GatingPreBool']
            args_balance = ['balance=Stimulus']
            args_fr = ['fr=WindowGatingClassify', 'fr=ConcatFactor2']
            args_counts_thr = ['counts_thr=15']
            args_area_list = ['area_list={0:s}'.format(area_list)]
            args_subject = ['subject={0:d}'.format(session)]
            args_area = ['area={0:s}'.format(area)]
            args_mode = ['mode=Normal']
            args_mode_seed = ['mode_seed=0']
            args_pseudo_seed = ['pseudo_seed=0']
            args_split = ['split=StratifiedStim']
            args_split_ind = ['split_ind=0']
            args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr,
                                                                 args_area_list, args_subject, args_area, args_mode,
                                                                 args_mode_seed, args_pseudo_seed, args_split, args_split_ind)))))

    for balance in ['StageGatingCenteredGatingOnly', 'StageGatingCenteredPostDist1Only']:
        for session in range(42):
            args_class = ['class=GatedStimulus']
            args_balance = ['balance={0:s}'.format(balance)]
            args_fr = ['fr=WindowMemoryClassify', 'fr=ConcatFactor2']
            args_counts_thr = ['counts_thr=15']
            args_area_list = ['area_list=PFC']
            args_subject = ['subject={0:d}'.format(session)]
            args_area = ['area=PFC']
            args_mode = ['mode=Normal']
            args_mode_seed = ['mode_seed=0']
            args_pseudo_seed = ['pseudo_seed=0']
            args_split = ['split=StratifiedBalanceSplit']
            args_split_ind = ['split_ind=0']
            args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr,
                                                                 args_area_list, args_subject, args_area, args_mode,
                                                                 args_mode_seed, args_pseudo_seed, args_split, args_split_ind)))))

    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))

    return args_version_from_job


main()






db = DataBase()