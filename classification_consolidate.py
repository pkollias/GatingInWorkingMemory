import sys
from rec_utils import *

def main():

    args_version = sys.argv[1:]

    # args_version = ['job_id=0', 'overwrite=True']

    # load analysis parameters
    version = job_scheduler(args_version, args_from_parse_func)

    classifier = ClassificationAnalysis(DataBase([]), version)
    db, md = classifier.db, classifier.db.md

    # overwrite check
    stem = classifier.get_train_test_stem()
    fname = 'cv_score_full'
    target_filename = classifier.get_path_base(fname, stem, cross=False)
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
    num_iters = num_shuffles

    # for every iteration load and save results in cv_score
    cv_score = np.zeros((43, 43, num_shuffles))
    for ii, shuffle in enumerate(shuffle_range):

        print(ii, num_iters)

        # class_list = version['class_list'].split('_')
        # balance_list = version['balance_list'].split('_')
        #
        # def version_modifier(class_i, balance_i):
        #     version_mod = version.copy()
        #     version_mod['class'] = class_i
        #     version_mod['balance'] = balance_i
        #     return version_mod

        score = {}
        # get training classifier

        shuffle_ind_i = shuffle_ind(shuffle, shuffle_range)
        fname = 'cv_score' if not shuffle else 'cv_score_{0:04d}'.format(shuffle)
        stem = classifier.get_train_test_stem()
        src_filename = classifier.get_path_base(fname, stem, cross=False)
        try:
            cv_score[:, :, shuffle_ind_i] = md.np_loader(src_filename)
        except:
            cv_score[:, :, shuffle_ind_i] = np.nan

    md.np_saver(cv_score, target_filename)


def args_from_parse_func(parse_version):

    args_version_list = []

    for split in ['WithinGroupTransfer']:
        for area_list, area in [('PFC', 'PFC'), ('Stri', 'Stri'), ('IT', 'IT')]:
            args_class_list = ['class=GatingPreBoolGeneralized']
            args_balance = ['balance=Stimulus']
            args_fr = ['fr=ConcatFactor2']
            args_counts_thr = ['counts_thr=12']
            args_area_list = ['area_list={0:s}'.format(area_list)]
            args_subject = ['subject=Gonzo_Oscar']
            args_area = ['area={0:s}'.format(area)]
            args_mode = ['mode=Normal']
            args_mode_seed = ['mode_seed=0']  # 100
            args_split = ['split={0:s}'.format(split)]
            args_pseudo_seed = ['pseudo_seed=0']
            args_shuffle_range = ['shuffle_range=1_250']
            args_version_list.extend(list(map(list, list(product(args_class_list, args_balance, args_fr, args_counts_thr,
                                                                 args_area_list, args_subject, args_area, args_mode,
                                                                 args_mode_seed, args_split, args_pseudo_seed, args_shuffle_range)))))

    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))

    return args_version_from_job


main()
