import sys
from rec_analyses import *
from rec import SignalSmoothing


def main():
    """ factor=, fr=, counts_thr=, area_list=, subject=, area=, mode=, mode_seed=, [overwrite=] """
    # args_version = ['factor=StimulusGatingPreBool', 'fr=ConcatFactor2', 'counts_thr=15',
    # 'area_list=PFC_Stri', 'subject=Gonzo_Oscar', 'area=PFC', 'mode=AreaShuffle', 'mode_seed=0']

    # load analysis parameters
    args = sys.argv
    args_version = args[1:]
    version = parse_vars(args_version)

    # create analysis object
    dpca = DemixedPrincipalComponent(DataBase([]), version)

    # overwrite check
    target_filename = [dpca.get_exec_filename(fn) for fn in ['fbt', 'dpca_obj', 'X_fit', 'X_tuple']]
    print(target_filename)
    if all([path.exists(fn) for fn in target_filename]) and ('overwrite' not in version.keys() or not version['overwrite']):
        exit()

    # load behavioral units from assemble
    pbt_full = dpca.db.md.np_loader(dpca.get_exec_filename('pbt'))

    # process
    # smoothen and crop data
    smoother = SignalSmoothing(signal.correlate, signal.windows.gaussian(9, 1.25))
    pbt = pbt_full.init_with_df(pbt_full.smooth_df(smoother))
    pbt.crop_timeseries(-50, 1000)
    dpca.pbt = pbt

    # exec
    dpca_obj, X_fit, X, X_trial, records, records_trial = dpca.exec(pbt)

    # convert to fbt
    fbt = {margin: FactorBehavioralTimeseries(fbt_df_from_dPCA(xf, records, dpca_obj, pbt.timebin_interval),
                                              pbt.condition_labels, pbt.timebin_interval)
           for margin, xf
           in X_fit.items()}

    for fn, var in zip(target_filename, [fbt, dpca_obj, X_fit, (X, X_trial, records, records_trial)]):
        dpca.db.md.np_saver(var, fn)


main()
