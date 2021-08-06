from rec import *


def main():

    # INIT
    md = MetaData()
    db = md.db_base_loader(['sessions', 'units', 'events', 'activity'])
    sessions, units, events, activity = db['sessions'], db['units'], db['events'], db['activity']

    # initialize units events dict
    units_events = {sess: pd.DataFrame() for sess in sessions.index}

    # for every session
    for sess in sessions.index:

        units_slice = units.loc[units['Session'].eq(sess)]
        events_slice = events.loc[events['Session'].eq(sess)]
        # initialize units_events dataframe for that session
        units_events_sess = pd.DataFrame(np.full((len(events_slice), len(units_slice)), False, dtype=bool),
                                         index=events_slice.index, columns=units_slice.index)

        # for every unit in session
        for unit_ind in units_slice.index:

            print(unit_ind)
            # load events data
            units_events_list = md.np_loader(md.preproc_dest_path(path.join('temp', 'units_events',
                                                                            'units_events_{0:s}_chan{1:03d}_unit{2:03d}.pkl'.format(*unit_ind))))
            # set events true
            for ue in units_events_list:
                units_events_sess.at[ue, unit_ind] = True

        # assign to units_events dict
        units_events[sess] = units_events_sess

    # SAVE
    md.np_saver(units_events, md.preproc_dest_path(path.join('units_events.pkl')))


main()
