import sys
from rec_analyses import *


def main():
    """ factor=, fr=, counts_thr=, mode=, area_list=, area=, subject=, mode_seed=, [overwrite=] """
    # args_version = ['factor=StimulusGating', 'fr=ConcatFactor2', 'counts_thr=20',
    # 'area_list=PFC_Stri', 'area=PFC', 'subject=Oscar', 'mode=AreaShuffle', 'mode_seed=0']

    # load analysis parameters
    args = sys.argv
    args_version = args[1:]
    version = parse_vars(args_version)

    # create analysis object
    dpca = DemixedPrincipalComponent(DataBase(['sessions', 'units', 'events', 'conditions']), version)

    # overwrite check
    target_filename = dpca.get_assemble_filename()
    print(target_filename)
    if path.exists(target_filename) and ('overwrite' not in version.keys() or not version['overwrite']):
        exit()

    units = dpca.db.tables['units']
    # create unit selection filters
    # valid units
    valid_units_events = dpca.db.md.np_loader(dpca.get_wrangle_filename())
    valid_units = valid_units_events.apply(bool)
    # single units
    single_units = units['UnitNum'].ne(0) & units['RatingCode'].ne(7)
    # area list units
    area_list_units = units['Area'].isin(dpca.version['area_list'].split('_'))
    # subject units
    sessions = dpca.db.tables['sessions']
    subject_sessions = sessions['Subject'].isin(dpca.version['subject'].split('_')).index
    subject_units = units['Session'].isin(subject_sessions)
    # apply filters
    dpca.db.tables['units'] = units.loc[valid_units_events & single_units & area_list_units & subject_units]

    dpca.db.md.np_saver(valid_units, target_filename)


main()
