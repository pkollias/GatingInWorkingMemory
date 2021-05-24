import sys
from rec_analyses import *


def main():
    """ factor=, fr=, counts_thr=, area_list=, subject=, [overwrite=] """
    # args_version = ['factor=StimulusGatingPreBool', 'fr=ConcatFactor2', 'counts_thr=15',
    # 'area_list=PFC_Stri', 'subject=Gonzo_Oscar']

    # load analysis parameters
    args = sys.argv
    args_version = args[1:]
    version = parse_vars(args_version)

    # create analysis object
    dpca = DemixedPrincipalComponent(DataBase(['sessions', 'units', 'events', 'conditions']), version)

    # overwrite check
    target_filename = dpca.get_filter_filename()
    print(target_filename)
    if path.exists(target_filename) and ('overwrite' not in version.keys() or not version['overwrite']):
        exit()

    # create unit selection filters
    units = dpca.db.tables['units']
    # valid units
    valid_units_events = dpca.db.md.np_loader(dpca.get_wrangle_filename())
    valid_units = valid_units_events.apply(bool)
    # single units
    single_units = units['UnitNum'].ne(0) & units['RatingCode'].ne(7)
    # area list units
    area_list_units = units['Area'].isin(dpca.version['area_list'].split('_'))
    # subject units
    sessions = dpca.db.tables['sessions']
    subject_sessions = sessions.loc[sessions['Subject'].isin(dpca.version['subject'].split('_'))].index
    subject_units = units['Session'].isin(subject_sessions)
    # apply filters
    dpca.db.tables['units'] = units.loc[valid_units & single_units & area_list_units & subject_units]

    # create behavioral_units table
    units_events = valid_units_events.loc[dpca.db.tables['units'].index]
    behavioral_units = pd.DataFrame([(sess, channum, unitnum, trialnum, stageindex)
                                     for (sess, channum, unitnum), events_list in units_events.iteritems()
                                     for _, trialnum, stageindex in events_list],
                                    columns=dpca.db.md.proc_imports['behavioral_units']['index'])

    dpca.db.md.np_saver(behavioral_units, target_filename)


main()
