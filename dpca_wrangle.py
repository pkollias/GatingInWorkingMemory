import sys
from rec_analyses import *


def main():


    # load analysis parameters
    args = sys.argv
    version_factor = args[1]
    version_fr = args[2]

    dpca = DemixedPrincipalComponent(DataBase(['units', 'events', 'conditions']), (version_factor, version_fr))

    target_filename = dpca.get_wrangle_filename()
    print(target_filename)
    if path.exists(target_filename):
        exit()

    units = dpca.db.tables['units']
    valid_units = units.apply(lambda row: bool(dpca.assess_unit_events(row.name)), axis=1)

    md.np_saver(valid_units, target_filename)


main()
