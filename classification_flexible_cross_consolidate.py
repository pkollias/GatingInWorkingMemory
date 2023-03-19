import sys
from rec_utils import *

def main():
    """ class_list, balance_list, fr, counts_thr, area_list, subject, area, mode, mode_seed, pseudo_seed, split,
    shuffle_range, [overwrite] """
    args_version = sys.argv[1:]
    # args_version = ['job_id=0', 'overwrite=True']
    version = job_scheduler(args_version, args_from_parse_func)

    classifier = ClassificationAnalysis(DataBase([]), version)
    db, md = classifier.db, classifier.db.md

    # overwrite check
    stem = (*classifier.get_train_test_stem(), 'consolidate')
    fname_stem = '_crossbin' if version['train_t'] == 'crossbin' else ''
    fname = 'cv_score_full' + fname_stem
    target_filename = classifier.get_path_base(fname, stem, cross=True)
    print(target_filename)
    if path.exists(target_filename) and ('overwrite' not in version.keys() or not eval(version['overwrite'])):
        exit()

    # util functions
    def range_str_to_range(range_str):
        range_lims = list(map(int, range_str.split('_')))
        return range(*range_lims)

    def shuffle_ind(shuffle_i, shuffle_range):
        return shuffle_i - shuffle_range[0]

    # get objects for which over to iterate for given version
    shuffle_range = range_str_to_range(version['shuffle_range'])
    num_shuffles = len(shuffle_range)
    classifier_ind_list = list(product(range(2), range(2)))
    split_split_ind_list = list(product(range(3), range(3)))
    subversions_iterator = list(product(shuffle_range, classifier_ind_list, split_split_ind_list))
    num_iters = len(subversions_iterator)

    # for every iteration load and save results in cv_score
    train_n = 1 if version['train_t'] == 'crossbin' else 43
    cv_score = np.zeros((4, 9, train_n, 43, num_shuffles))
    for ii, (shuffle, classifier_ind, split_split_ind) in enumerate(subversions_iterator):

        print(ii, num_iters)
        class_ravel_ind = np.ravel_multi_index(classifier_ind, (2, 2))
        split_ravel_ind = np.ravel_multi_index(split_split_ind, (3, 3))

        class_list = version['class_list'].split('_')
        balance_list = version['balance_list'].split('_')

        def version_modifier(class_i, balance_i):
            version_mod = version.copy()
            version_mod['class'] = class_i
            version_mod['balance'] = balance_i
            return version_mod

        # create analysis object
        classifier_list = [ClassificationAnalysis(DataBase([]), version_modifier(class_i, balance_i)) for
                           class_i, balance_i
                           in zip(class_list, balance_list)]
        score = {}
        # get training classifier
        class_train = classifier_list[classifier_ind[0]]
        db_train, md_train = class_train.db, class_train.db.md
        train_dict = md_train.np_loader(class_train.get_path_base('train_test', class_train.get_train_test_stem()))
        # get split of splits
        if class_train.version['split'] == 'StratifiedBalanceSplit_StimHalf':
            train_key = list(np.unique(list(map(itemgetter(0), train_dict.keys()))))[split_split_ind[0]]
        else:
            train_key = list(train_dict.keys())[split_split_ind[0]]

        # get testing classifier
        class_test = classifier_list[classifier_ind[1]]
        db_test, md_test = class_test.db, class_test.db.md
        test_dict = md_test.np_loader(class_test.get_path_base('train_test', class_test.get_train_test_stem()))
        # get split of splits
        if class_test.version['split'] == 'StratifiedBalanceSplit_StimHalf':
            test_key = list(np.unique(list(map(itemgetter(0), test_dict.keys()))))[split_split_ind[1]]
        else:
            test_key = list(test_dict.keys())[split_split_ind[1]]

        stem = (*classifier_list[0].get_train_test_stem(),
                'intermediate',
                '_'.join([class_train.version['class'], train_key, class_test.version['class'], test_key]))

        shuffle_ind_i = shuffle_ind(shuffle, shuffle_range)
        fname_stem = '_crossbin' if version['train_t'] == 'crossbin' else ''
        fname = ('cv_score' if not shuffle_ind_i else 'cv_score_{0:04d}'.format(shuffle_ind_i)) + fname_stem
        src_filename = classifier_list[0].get_path_base(fname, stem, cross=True)
        try:
            cv_score[class_ravel_ind, split_ravel_ind, :, :, shuffle_ind_i] = md.np_loader(src_filename)
        except:
            cv_score[class_ravel_ind, split_ravel_ind, :, :, shuffle_ind_i] = np.nan

    md.np_saver(cv_score, target_filename)


def args_from_parse_func(parse_version):

    args_version_list = []

    for split, shuffle_range in [('StratifiedBalanceSplit_StimHalf', '1_101')]: #351
        for area_list, area in [('PFC', 'PFC')]: #Stri, IT
            args_class_list = ['class_list=Stimulus_GatedStimulus']
            args_balance = ['balance_list=StageGatingPrePostMemory_StageGatingPrePostSensory']
            args_fr = ['fr=ConcatFactor2']
            args_counts_thr = ['counts_thr=12']
            args_area_list = ['area_list={0:s}'.format(area_list)]
            args_subject = ['subject=Gonzo_Oscar']
            args_area = ['area={0:s}'.format(area)]
            args_mode = ['mode=Normal']
            args_mode_seed = ['mode_seed=0']
            args_split = ['split=StratifiedBalanceSplit_StimHalf']
            args_pseudo_seed = ['pseudo_seed=0']
            args_train_t = ['train_t=crossbin']
            args_shuffle_range = ['shuffle_range={0:s}'.format(shuffle_range)]
            args_version_list.extend(list(map(list, list(product(args_class_list, args_balance, args_fr, args_counts_thr,
                                                                 args_area_list, args_subject, args_area, args_mode,
                                                                 args_mode_seed, args_split, args_pseudo_seed, args_train_t,
                                                                 args_shuffle_range)))))

    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))

    return args_version_from_job


main()
