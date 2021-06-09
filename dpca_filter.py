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
    db, md = dpca.db, dpca.db.md

    # overwrite check
    target_filename = dpca.get_path_base('filter', dpca.get_filter_stem())
    print(target_filename)
    if path.exists(target_filename) and ('overwrite' not in version.keys() or not eval(version['overwrite'])):
        exit()

    # create unit selection filters
    units = db.tables['units']
    # valid units
    valid_units_events = md.np_loader(dpca.get_path_base('valid_units', dpca.get_wrangle_stem()))
    valid_units = valid_units_events.apply(bool)
    # single units
    single_units = units['UnitNum'].ne(0) & units['RatingCode'].ne(7)
    # area list units
    area_list_units = units['Area'].isin(dpca.version['area_list'].split('_'))
    # subject units
    sessions = db.tables['sessions']
    subject_sessions = sessions.loc[sessions['Subject'].isin(dpca.version['subject'].split('_'))].index
    subject_units = units['Session'].isin(subject_sessions)
    # apply filters
    db.tables['units'] = units.loc[valid_units & single_units & area_list_units & subject_units]

    # create behavioral_units table
    units_events = valid_units_events.loc[db.tables['units'].index]
    behavioral_units = pd.DataFrame([(sess, channum, unitnum, trialnum, stageindex)
                                     for (sess, channum, unitnum), events_list in units_events.iteritems()
                                     for _, trialnum, stageindex in events_list],
                                    columns=md.proc_imports['behavioral_units']['index'])

    md.np_saver(behavioral_units, target_filename)


main()
