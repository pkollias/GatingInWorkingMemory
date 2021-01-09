from metadata import *

def trials_in_unit(trials, activity, unit_ind):

    sess = unit_ind[0]
    activity_collection = activity.loc[unit_ind]

    return trials_in_segment(trials, sess, activity_collection)



def trials_in_segment(trials, sess, activity_collection):

    trials_slice = trials.loc[sess]

    row_in_segment = lambda row: timestamp_interval_within_activity(row['StartTimeStamp'], row['EndTimeStamp'], activity_collection)
    trials_slice_filter = trials_slice.apply(row_in_segment, axis=1)
    trialnum_list = trials_slice[trials_slice_filter].index

    return [(sess, trialnum) for trialnum in trialnum_list]



def timestamp_interval_within_activity(interval_start, interval_end, activity_collection):

    return (activity_collection['SegmentStart'].le(interval_start) &
            activity_collection['SegmentEnd'].gt(interval_start)).any() \
           and \
           (activity_collection['SegmentStart'].le(interval_end) &
            activity_collection['SegmentEnd'].gt(interval_end)).any()
