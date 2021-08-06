import sys
from rec import *
from rec_db import timestamp_interval_within_activity


def main():

    # parameters
    args = sys.argv
    u_iloc = int(args[1])

    # INIT
    md = MetaData()
    db = md.db_base_loader(['units', 'events', 'activity'])
    units, events, activity = db['units'], db['events'], db['activity']

    sm = SamplingMethod()

    unit_ind = units.iloc[u_iloc].name
    sess, _, _ = unit_ind
    events_slice = events.loc[events['Session'].eq(sess)]
    activity_collection = activity.loc[unit_ind]

    def event_within_activity(event_row):
        if np.isnan(event_row['StageDuration']):
            return False
        event_onset_point = SamplingPoint(sm, stamp=event_row['StageTimeIndex'])
        event_window = SamplingInterval(SamplingPoint(sm, t=0 * qu.ms),
                                        SamplingPoint(sm, t=event_row['StageDuration'] * qu.ms))
        event_fr_interval = SamplingInterval(event_onset_point, event_onset_point.copy()).get_offset(event_window)
        return timestamp_interval_within_activity(event_fr_interval.start.stamp, event_fr_interval.end.stamp,
                                                  activity_collection)

    # for every event
    unit_mask = events_slice.apply(event_within_activity, axis=1)
    unit_events_list = unit_mask.loc[unit_mask].index

    # SAVE
    md.np_saver(unit_events_list,
                md.preproc_dest_path(path.join('temp', 'units_events',
                                               'units_events_{0:s}_chan{1:03d}_unit{2:03d}.pkl'.format(*unit_ind))))


main()
