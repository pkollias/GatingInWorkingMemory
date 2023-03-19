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

    # util functions
    def range_str_to_range(range_str):
        range_lims = list(map(int, range_str.split('_')))
        return range(*range_lims)

    # get objects for which over to iterate for given version
    mode_seed_range = range_str_to_range(version['mode_seed_range'])
    num_mode_seeds = len(mode_seed_range)

    cv_score = np.zeros((73, num_mode_seeds))
    for mode_seed_i, mode_seed in enumerate(mode_seed_range):

        try:
            version['mode_seed'] = '{0:03d}'.format(mode_seed)

            classifier_ind = list(map(int, version['classifier_ind'].split('_')))
            split_split_ind = list(map(lambda x: [int(el) for el in list(x)], version['split_split_ind'].split('_')))
            class_list = version['class_list'].split('_')
            balance_list = version['balance_list'].split('_')

            def version_modifier(class_i, balance_i):
                version_mod = version.copy()
                version_mod['class'] = class_i
                version_mod['balance'] = balance_i
                return version_mod

            # create analysis object
            classifier_list = [ClassificationAnalysis(DataBase([]), version_modifier(class_i, balance_i))
                               for class_i, balance_i
                               in zip(class_list, balance_list)]

            # get training classifier
            class_train = classifier_list[classifier_ind[0]]
            train_dict = md.np_loader(class_train.get_path_base('train_test', class_train.get_train_test_stem()))
            # get split of splits
            train_key_abstract = itemgetter(*split_split_ind[0])(
                list(np.unique(list(map(itemgetter(0), train_dict.keys())))))  ### 2GROUP
            train_key_abstract = list(train_key_abstract) if type(train_key_abstract) is tuple else [train_key_abstract]

            # get testing classifier
            class_test = classifier_list[classifier_ind[1]]
            test_dict = md.np_loader(class_test.get_path_base('train_test', class_test.get_train_test_stem()))
            # get split of splits
            test_key_abstract = itemgetter(*split_split_ind[1])(
                list(np.unique(list(map(itemgetter(0), test_dict.keys())))))  ### 2GROUP
            test_key_abstract = list(test_key_abstract) if type(test_key_abstract) is tuple else [test_key_abstract]

            # overwrite check
            stem = (*classifier_list[0].get_train_test_stem(),
                    'intermediate',
                    '_'.join(
                        [class_train.version['class'], ''.join(sorted(train_key_abstract)), class_test.version['class'],
                         ''.join(sorted(test_key_abstract))]))  ### 2GROUP
            fname_stem = '_{0:s}'.format(version['train_t']) if version['train_t'] in ['crossbin', 'sens', 'mem',
                                                                                       'diagonal'] else ''
            fname = ('cv_score' if not int(version['shuffle']) else 'cv_score_{0:04d}'.format(
                int(version['shuffle']))) + fname_stem

            src_filename = classifier_list[0].get_path_base(fname, stem, cross=True)

            cv_score[:, mode_seed_i] = md.np_loader(src_filename).squeeze()
        except:
            cv_score[:, mode_seed_i] = np.nan

    # overwrite check ### PK This should normally be before script run
    del version['mode_seed'] ### PK Still Hacking
    version['mode_seed'] = version['mode_seed_range']

    stem = (*classifier.get_train_test_stem(),
            'mode_average',
            '_'.join([class_train.version['class'],
                      ''.join(sorted(train_key_abstract)),
                      class_test.version['class'],
                      ''.join(sorted(test_key_abstract))]))  ### 2GROUP
    fname_stem = version['train_t']
    fname = 'cv_score_mean_' + fname_stem
    target_filename = classifier.get_path_base(fname, stem, cross=True)
    print(target_filename)
    if path.exists(target_filename) and ('overwrite' not in version.keys() or not eval(version['overwrite'])):
        exit()

    md.np_saver(cv_score, target_filename)


def args_from_parse_func(parse_version):

    args_version_list = []

    for (classifier_ind_arg, test_split_arg) in product([('0', '0')], ['2', '0', '1']):
        for (split_split_ind_arg, train_t) in [((test_split_arg, test_split_arg), 'diagonal'),
                                               (('012', test_split_arg), 'sens'),
                                               (('12', test_split_arg), 'mem'),
                                               (('0', test_split_arg), 'mem')]:
            for mode_seed_range in [('0_1000')]:
                for area_list, area in [('PFC', 'PFC')]:
                    args_class_list = ['class_list=Stimulus_GatedStimulus']
                    args_balance = ['balance_list=StageGatingPrePostMemory_StageGatingPrePostSensory']
                    args_fr = ['fr=ConcatFactor2']
                    args_counts_thr = ['counts_thr=12']
                    args_area_list = ['area_list={0:s}'.format(area_list)]
                    args_subject = ['subject=Gonzo_Oscar']
                    args_area = ['area={0:s}'.format(area)]
                    args_mode = ['mode=Normal']
                    args_mode_seed_range = ['mode_seed_range={0:s}'.format(mode_seed_range)]
                    args_split = ['split=StratifiedBalanceSplit_StimHalf']
                    args_pseudo_seed = ['pseudo_seed=0']
                    args_train_t = ['train_t={0:s}'.format(train_t)]
                    args_shuffle = ['shuffle=0']
                    args_classifier_ind = ['classifier_ind={0:s}_{1:s}'.format(*classifier_ind_arg)]
                    args_split_split_ind = ['split_split_ind={0:s}_{1:s}'.format(*split_split_ind_arg)]
                    args_version_list.extend(
                        list(map(list, list(product(args_class_list, args_balance, args_fr, args_counts_thr,
                                                    args_area_list, args_subject, args_area, args_mode,
                                                    args_mode_seed_range, args_split, args_pseudo_seed, args_train_t,
                                                    args_shuffle, args_classifier_ind, args_split_split_ind)))))

    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))

    return args_version_from_job


main()
