import sys
from rec_analyses import *


def main():
    """ factor=, fr=, counts_thr=, area_list=, subject=, area=, mode=, mode_seed=, [overwrite=] """
    # args_version = ['factor=StimulusGatingPreBool', 'fr=ConcatFactor2', 'counts_thr=20',
    # 'area_list=PFC_Stri', 'subject=Gonzo_Oscar', 'area=PFC', 'mode=AreaShuffle', 'mode_seed=0']

    # load analysis parameters
    args = sys.argv
    args_version = args[1:]
    version = parse_vars(args_version)

    # create analysis object
    dpca = DemixedPrincipalComponent(DataBase([]), version)
    db, md = dpca.db, dpca.db.md

    # overwrite check
    target_filename = dpca.get_path_base('significance', dpca.get_exec_stem())
    print(target_filename)
    if path.exists(target_filename) and ('overwrite' not in version.keys() or not eval(version['overwrite'])):
        exit()

    src_filename = [dpca.get_path_base(fn, dpca.get_exec_stem()) for fn in ['dpca_obj', 'X_tuple']]
    dpca_obj, (X, X_trial, _, _) = tuple(
        [md.np_loader(fn) for fn in src_filename])

    # process
    # run cross-validated mean classification score significance analyses
    significance = dpca_obj.significance_analysis(X, X_trial, n_shuffles=100, n_splits=20, axis=True)

    md.np_saver(significance, target_filename)


main()
