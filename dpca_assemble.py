import sys
from rec_utils import *


def main():
    """ factor=, fr=, counts_thr=, area_list=, subject=, area=, mode=, mode_seed=, [overwrite=] """
    # args_version = ['factor=StimulusGatingPreBool', 'fr=ConcatFactor2', 'counts_thr=15',
    # 'area_list=PFC_Stri', 'subject=Gonzo_Oscar', 'area=PFC', 'mode=AreaShuffle', 'mode_seed=0']

    # load analysis parameters
    args = sys.argv
    args_version = args[1:]
    version = parse_vars(args_version)

    # create analysis object
    dpca = DemixedPrincipalComponent(DataBase(['units', 'events', 'conditions']), version)
    db, md = dpca.db, dpca.db.md

    # overwrite check
    target_filename = dpca.get_path_base('pbt', dpca.get_exec_stem())
    print(target_filename)
    if path.exists(target_filename) and ('overwrite' not in version.keys() or not eval(version['overwrite'])):
        exit()

    # load behavioral units from filter
    behavioral_units_filter = md.np_loader(dpca.get_path_base('filter', dpca.get_filter_stem()))

    # init params
    units_index = md.preproc_imports['units']['index']
    events_index = md.preproc_imports['events']['index']
    mode_seed = int(dpca.version['mode_seed'])
    area = dpca.version['area']
    v_factor_params = factor_generate_conditions(dpca.version['factor'])
    condition_columns = v_factor_params['condition_columns']
    counts_thr = int(dpca.version['counts_thr'])

    # process
    # get behavioral unit indices
    buf_units = list(set(zip_columns(behavioral_units_filter, units_index)))
    # filter units table with behavioral units
    db.tables['units'] = db.tables['units'].loc[buf_units]
    units = db.tables['units']

    if dpca.version['mode'] == 'AreaShuffle':
        # shuffle area labels and take shuffled area units indices
        units['AreaShuffle'] = list(units['Area'].sample(frac=1, random_state=mode_seed)) if mode_seed else units['Area']
        area_units = units.loc[units['AreaShuffle'].eq(area)].index
    elif dpca.version['mode'] == 'Bootstrap':
        # get min size of areas, calculate fraction to subsample, get area units indices subsampled (constant seed)
        area_counts = units['Area'].value_counts()
        min_size = area_counts.loc[list_param_to_list(dpca.version['area_list'])].min()
        area_units = units.loc[units['Area'].eq(area)].sample(min_size, random_state=0).index
    elif dpca.version['mode'] == 'Full':
        area_units = units.loc[units['Area'].eq(area)].index

    # TODO: modifications for replacement sampling of units, correct later
    area_units_df = pd.DataFrame(area_units, columns=['Unit'])
    area_units_df = unzip_columns(area_units_df, 'Unit', units_index)
    area_units_df['Unit_Code'] = area_units_df.index

    # isolate behavioral units of area selection
    area_behavioral_units = pd.merge(behavioral_units_filter, area_units_df, on=units_index)
    # merge with condition columns
    events_conditions = db.tables['events_conditions'][condition_columns].reset_index()
    area_behavioral_units_conditions = pd.merge(area_behavioral_units, events_conditions, on=events_index)
    # groupby units x conditions
    area_conditions_grouper = area_behavioral_units_conditions.groupby(units_index + ['Unit_Code'] + condition_columns)
    # sample counts_thr events for every unit x condition (replace for Bootstrap)
    behavioral_units = area_conditions_grouper.sample(counts_thr,
                                                      random_state=mode_seed,
                                                      replace=dpca.version['mode'] == 'Bootstrap')

    pbt = pbt_from_behavioral_units(condition_columns, version['fr'], behavioral_units, db)

    md.np_saver(pbt, target_filename)


main()
