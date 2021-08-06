import sys
from rec_utils import *
from rec import SignalSmoothing


def main():

    args_version = sys.argv[1:]

    """ factor=, fr=, counts_thr=, area_list=, subject=, area=, mode=, mode_seed=, [overwrite=] """
    # args_version = ['factor=StimulusGatingPreBool', 'fr=ConcatFactor2', 'counts_thr=20',
    #                 'area_list=PFC_Stri', 'subject=Gonzo_Oscar', 'area=PFC', 'mode=AreaShuffle', 'mode_seed=0']
    # args_version = ['job_id=0', 'overwrite=True']

    # load analysis parameters
    version = job_scheduler(args_version, args_from_parse_func)

    # create analysis object
    dpca = DemixedPrincipalComponent(DataBase([]), version)
    db, md = dpca.db, dpca.db.md

    # overwrite check
    target_filename = [dpca.get_path_base(fn, dpca.get_exec_stem()) for fn in ['fbt', 'dpca_obj', 'X_fit', 'X_tuple']]
    print(target_filename)
    if all([path.exists(fn) for fn in target_filename]) and ('overwrite' not in version.keys() or not eval(version['overwrite'])):
        exit()

    # load behavioral units from assemble
    pbt_full = md.np_loader(dpca.get_path_base('pbt', dpca.get_exec_stem()))

    # process
    # smoothen and crop data
    smoother = SignalSmoothing(signal.correlate, signal.windows.gaussian(9, 1.25))
    pbt = pbt_full.init_with_df(pbt_full.smooth_df(smoother))
    pbt.crop_timeseries(-50, 1000)
    dpca.pbt = pbt

    # exec
    dpca_obj, X_fit, X, X_trial, records, records_trial = dpca.exec()

    # convert to fbt
    fbt = {margin: FactorBehavioralTimeseries(fbt_df_from_dPCA(xf, records, dpca_obj, pbt.timebin_interval),
                                              pbt.condition_labels, pbt.timebin_interval)
           for margin, xf
           in X_fit.items()}

    for fn, var in zip(target_filename, [fbt, dpca_obj, X_fit, (X, X_trial, records, records_trial)]):
        md.np_saver(var, fn)
    # fbt, dpca_obj, X_fit, (X, X_trial, records, records_trial) = [md.np_loader(fn) for fn in target_filename]


def args_from_parse_func(parse_version):

    args_version_list = []

    for area_list, area in [('PFC', 'PFC'), ('Stri', 'Stri'), ('IT', 'IT')]:
        for subject in ['Gonzo', 'Oscar', 'Gonzo_Oscar']:
            for factor, counts_thr in [('GatedStimulusPostDistMemory', '6'), ('GatedStimulusPostDistMemory', '8'),
                                       ('GatedStimulusPostDistMemory', '10'), ('GatedStimulusPostDistMemory', '12')]:
                args_factor = ['factor={0:s}'.format(factor)]
                args_fr = ['fr=ConcatFactor2']
                args_counts_thr = ['counts_thr={0:s}'.format(counts_thr)]
                args_area_list = ['area_list={0:s}'.format(area_list)]
                args_subject = ['subject={0:s}'.format(subject)]
                args_area = ['area={0:s}'.format(area)]
                args_mode = ['mode=Full']
                args_mode_seed = ['mode_seed={0:d}'.format(ii) for ii in range(1)]
                args_version_list.extend(list(map(list, list(product(args_factor, args_fr, args_counts_thr, args_area_list,
                                                                     args_subject, args_area, args_mode, args_mode_seed)))))

    for area_list, area in [('PFC', 'PFC'), ('Stri', 'Stri'), ('IT', 'IT')]:
        for subject in ['Gonzo', 'Oscar', 'Gonzo_Oscar']:
            for counts_thr in ['10', '12', '16']:
                for factor in ['GatPostStimulusRuleStim', 'TargPostStimulusRuleStim',
                               'GatingPreBool', 'StimulusGating', 'StimulusGatingPreBool']:
                    args_factor = ['factor={0:s}'.format(factor)]
                    args_fr = ['fr=ConcatFactor2']
                    args_counts_thr = ['counts_thr={0:s}'.format(counts_thr)]
                    args_area_list = ['area_list={0:s}'.format(area_list)]
                    args_subject = ['subject={0:s}'.format(subject)]
                    args_area = ['area={0:s}'.format(area)]
                    args_mode = ['mode=Full']
                    args_mode_seed = ['mode_seed={0:d}'.format(ii) for ii in range(1)]
                    args_version_list.extend(list(map(list, list(product(args_factor, args_fr, args_counts_thr, args_area_list,
                                                                         args_subject, args_area, args_mode, args_mode_seed)))))

    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))

    return args_version_from_job


main()
