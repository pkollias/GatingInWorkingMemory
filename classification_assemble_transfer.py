import sys
from rec_utils import *


def main():

    args_version = sys.argv[1:]

    """ class=, balance, fr=, counts_thr=, area_list=, subject=, area=, mode=, mode_seed=, [overwrite=] """
    # args_version = ['class=Stimulus', 'balance=StageGatingCentered', 'fr=ConcatFactor2', 'counts_thr=30',
    # 'area_list=PFC_Stri', 'subject=Gonzo_Oscar', 'area=PFC', 'mode=Normal', 'mode_seed=0']
    # args_version = ['job_id=0']

    # load analysis parameters
    version = job_scheduler(args_version, args_from_parse_func)

    # create analysis object
    classifier = ClassificationAnalysis(DataBase(['units', 'events', 'conditions']), version)
    db, md = classifier.db, classifier.db.md

    # overwrite check
    target_filename = classifier.get_path_base('pbt', classifier.get_assemble_stem())
    print(target_filename)
    if path.exists(target_filename) and ('overwrite' not in version.keys() or not eval(version['overwrite'])):
        exit()

    # create unit selection filters
    behavioral_units_filter = md.np_loader(classifier.get_path_base('filter', classifier.get_filter_stem()))

    # init params
    units_index = md.preproc_imports['units']['index']
    events_index = md.preproc_imports['events']['index']
    mode_seed = int(classifier.version['mode_seed'])
    condition_columns = classifier.get_assembly_condition_columns()

    # process
    # get behavioral unit indices

    # filter units table with behavioral units
    buf_units = zip_columns(behavioral_units_filter, units_index).unique()
    db.tables['units'] = db.tables['units'].loc[buf_units]

    # transfer area units from transfer analysis
    version_transfer = version.copy()
    version_transfer['class'] = version['class_transfer']
    version_transfer['balance'] = version['balance_transfer']
    version_transfer['counts_thr'] = version['counts_thr_transfer']
    classifier_transfer = ClassificationAnalysis(DataBase(['units', 'events', 'conditions']), version)
    pbt_transfer = md.np_loader(classifier_transfer.get_path_base('pbt', classifier_transfer.get_assemble_stem()))
    pbt_transfer_area_units_df = pbt_transfer.df[['Unit', 'Unit_Code']].drop_duplicates()
    unzipped_unit = unzip_columns(pbt_transfer_area_units_df, 'Unit', md.preproc_imports['units']['index'])
    pbt_transfer_area_units_unzipped_df = pd.concat([unzipped_unit, pbt_transfer_area_units_df['Unit_Code']], axis=1)
    area_units_df = pbt_transfer_area_units_unzipped_df.sort_values('Unit_Code').reset_index(drop=True)


    # isolate behavioral units of area selection
    area_behavioral_units = pd.merge(behavioral_units_filter, area_units_df, on=units_index)
    # merge with condition columns
    events_conditions = db.tables['events_conditions'][condition_columns].reset_index()
    area_behavioral_units_conditions = pd.merge(area_behavioral_units, events_conditions, on=events_index)
    # groupby units x conditions
    area_conditions_grouper = area_behavioral_units_conditions.groupby(units_index + ['Unit_Code'] + condition_columns)
    # sample counts_thr events for every unit x condition (replace for Bootstrap)
    replace_events = classifier.version['mode'] in ['Bootstrap', 'BootstrapEvents']
    seed_events = mode_seed
    counts_thr = area_conditions_grouper.apply(len).min()
    behavioral_units = area_conditions_grouper.sample(counts_thr,
                                                      random_state=get_seed(('events', seed_events)),
                                                      replace=replace_events)

    pbt = pbt_from_behavioral_units(condition_columns, version['fr'], behavioral_units, db)

    md.np_saver(pbt, target_filename)


def args_from_parse_func(parse_version):

    args_version_list = []

    for area_list, area in [('PFC', 'PFC'), ('Stri', 'Stri'), ('IT', 'IT')]:
        args_class = ['class=GatingPreBoolGeneralizedCue']
        args_balance = ['balance=Cue']
        args_fr = ['fr=ConcatFactor2']
        args_counts_thr = ['counts_thr={0:s}'.format(counts_thr) for counts_thr in ['1']]
        args_area_list = ['area_list={0:s}'.format(area_list)]
        args_subject = ['subject=Gonzo_Oscar']
        args_area = ['area={0:s}'.format(area)]
        args_mode = ['mode=Normal']
        args_mode_seed = ['mode_seed={0:s}'.format(mode_seed) for mode_seed in [str(ms) for ms in range(1)]]
        args_class_transfer = ['class_transfer=GatingPreBoolGeneralized']
        args_balance_transfer = ['balance_transfer=Stimulus']
        args_counts_thr_transfer = ['counts_thr_transfer=12']
        args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr,
                                                             args_area_list, args_subject, args_area,
                                                             args_mode, args_mode_seed,
                                                             args_class_transfer, args_balance_transfer, args_counts_thr_transfer)))))

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
    #         args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr,
    #                                                              args_area_list, args_subject, args_area,
    #                                                              args_mode, args_mode_seed)))))

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
    #         args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr,
    #                                                              args_area_list, args_subject, args_area,
    #                                                              args_mode, args_mode_seed)))))
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
    #         args_version_list.extend(list(map(list, list(product(args_class, args_balance, args_fr, args_counts_thr,
    #                                                              args_area_list, args_subject, args_area,
    #                                                              args_mode, args_mode_seed)))))

    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))

    return args_version_from_job



main()
