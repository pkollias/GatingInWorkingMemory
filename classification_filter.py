import sys
from rec_analyses import *


def main():
    """ class=, balance, fr=, counts_thr=, area_list=, subject=, [overwrite=] """
    # args_version = ['class=GatingPreBool', 'balance=Stimulus', 'fr=ConcatFactor2', 'counts_thr=15',
    # 'area_list=PFC_Stri', 'subject=Gonzo_Oscar']

    # load analysis parameters
    args = sys.argv
    args_version = args[1:]
    version = parse_vars(args_version)

    # create analysis object
    classifier = ClassificationAnalysis(DataBase(['sessions', 'units', 'events', 'conditions']), version)

    # overwrite check
    target_filename = classifier.get_filter_filename()
    print(target_filename)
    if path.exists(target_filename) and ('overwrite' not in version.keys() or not version['overwrite']):
        exit()

    # create unit selection filters
    units = classifier.db.tables['units']
    # valid units
    valid_units_events = classifier.db.md.np_loader(classifier.get_wrangle_filename())
    valid_units = valid_units_events.apply(bool)
    # single units
    single_units = units['UnitNum'].ne(0) & units['RatingCode'].ne(7)
    # area list units
    area_list = classifier.version['area_list'].split('_')
    area_list_units = units['Area'].isin(area_list)
    # subject units
    sessions = classifier.db.tables['sessions']
    subject_sessions = sessions.loc[sessions['Subject'].isin(classifier.version['subject'].split('_'))].index
    subject_units = units['Session'].isin(subject_sessions)
    # apply filters
    classifier.db.tables['units'] = units.loc[valid_units & single_units & area_list_units & subject_units]

    # create behavioral_units table
    units_events = valid_units_events.loc[classifier.db.tables['units'].index]
    behavioral_units = pd.DataFrame([(sess, channum, unitnum, trialnum, stageindex)
                                     for (sess, channum, unitnum), events_list in units_events.iteritems()
                                     for _, trialnum, stageindex in events_list],
                                    columns=classifier.db.md.proc_imports['behavioral_units']['index'])

    classifier.db.md.np_saver(behavioral_units, target_filename)


main()
