import sys
from rec_analyses import *


def main():
    """ factor=, fr=, counts_thr=, fr_thr=, [overwrite=] """
    # args_version = ['factor=StimulusGatingPreBool', 'fr=ConcatFactor2', 'counts_thr=15', 'fr_thr=100']

    # load analysis parameters
    args = sys.argv
    args_version = args[1:]
    version = parse_vars(args_version)

    # create analysis object
    dpca = DemixedPrincipalComponent(DataBase(['trials', 'units', 'events', 'conditions']), version)
    db, md = dpca.db, dpca.db.md

    # overwrite check
    target_filename = dpca.get_path_base('valid_units', dpca.get_wrangle_stem())
    print(target_filename)
    if path.exists(target_filename) and ('overwrite' not in version.keys() or not eval(version['overwrite'])):
        exit()

    # mark units as valid based on number of events and
    units = db.tables['units']
    valid_units_events = units.apply(lambda row: dpca.assess_unit_events(row.name), axis=1)

    md.np_saver(valid_units_events, target_filename)


main()
