import sys
from rec_analyses import *


def main():
    """ class=, balance=, fr=, counts_thr=, [overwrite=] """
    # args_version = ['class=GatingPreBool', 'balance=Stimulus', 'fr=ConcatFactor2', 'counts_thr=15']

    # load analysis parameters
    args = sys.argv
    args_version = args[1:]
    version = parse_vars(args_version)

    # create analysis object
    classifier = ClassificationAnalysis(DataBase(['trials', 'units', 'events', 'conditions']), version)

    # overwrite check
    target_filename = classifier.get_wrangle_filename()
    print(target_filename)
    if path.exists(target_filename) and ('overwrite' not in version.keys() or not version['overwrite']):
        exit()

    # mark units as valid based on number of events and
    classifier.db.merge_events_conditions_trials()
    units = classifier.db.tables['units']
    valid_units_events = units.apply(lambda row: classifier.assess_unit_events(row.name), axis=1)

    classifier.db.md.np_saver(valid_units_events, target_filename)


main()
