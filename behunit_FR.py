import sys
from rec import *
from aov_format import *


def main():

    # parameters
    args = sys.argv
    u_iloc = int(args[1])
    version_fr = args[2]
    v_fr_params = anova_version_fr_params(version_fr)
    t_start = v_fr_params['t_start']
    t_end = v_fr_params['t_end']
    timebin = v_fr_params['timebin']
    timestep = v_fr_params['timestep']
    event_mask = v_fr_params['event_mask']

    augment_columns = ['StageTimeIndex']

    # init tables and slices
    md = MetaData()
    db = md.db_base_loader(['sessions', 'trials', 'events', 'units', 'conditions', 'activity'])
    sessions, trials, events, units, conditions, activity =\
        db['sessions'], db['trials'], db['events'], db['units'], db['conditions'], db['activity']
    unit_entry = units.iloc[u_iloc]
    unit_ind = tuple(unit_entry[md.preproc_imports['units']['index']])
    sess, channum, unitnum = unit_ind
    target_filename = md.proc_dest_path(path.join('BehavioralUnits', 'FiringRates',
                                                  behunit_params_str(version_fr, timebin, timestep, t_start, t_end)),
                                        'behunit_FR_{0:s}_chan{1:03d}_unit{2:03d}.pkl'.format(sess, channum, unitnum))
    print(target_filename)
    if path.exists(target_filename):
        exit()

    trials_index = md.preproc_imports['trials']['index']
    events_index = md.preproc_imports['events']['index']

    # filter trials
    trials_slice = trials.loc[sess].reset_index(drop=True)

    # filter firing rate version only events
    events_slice = events.loc[sess].reset_index(drop=True)
    conditions_slice = conditions.loc[conditions['Session'].eq(sess)]
    events_conditions = pd.merge(events_slice, conditions_slice, on=events_index)
    events_conditions_slice = filter_df_wrapper(events_conditions, event_mask['column'], event_mask['wrapper'], event_mask['arg'])['df']

    # merge into final df
    trials_columns = trials_index + list(trials_slice.columns.difference(events_conditions_slice.columns))
    df = pd.merge(events_conditions_slice, trials_slice[trials_columns], on=trials_index)
    df.drop(df.columns.difference(events_index + augment_columns), axis=1, inplace=True)
    df.set_index(events_index, inplace=True, drop=False)
    df.sort_index(inplace=True)

    # preprocessing
    sm = SamplingMethod()
    event_window = SamplingInterval(SamplingPoint(sm, t=t_start * qu.ms),
                                    SamplingPoint(sm, t=t_end * qu.ms))

    # create events list for every timebin level
    timebin_fr_dict = {}

    session = Session(sessions.loc[sess], sm)
    data_session = md.np_loader(md.spike_dest_path('ts', sess, channum, unitnum))
    session_spiketrain = SpikeTrain(sess, sm, data_session, session.session_interval, src_data=data_session)
    activity_collection = activity.loc[unit_ind]

    # for every event
    for event_entry in df.itertuples():

        event_onset_point = SamplingPoint(sm, stamp=event_entry.StageTimeIndex)
        event_fr_interval = SamplingInterval(event_onset_point, event_onset_point.copy()).get_offset(event_window)

        # if event is within activity window
        if timestamp_interval_within_activity(event_fr_interval.start.stamp, event_fr_interval.end.stamp, activity_collection):

            # estimate firing rate for every timebin
            slice_fr_interval = session_spiketrain.index_of_slice(event_fr_interval)
            slice_fr_time_bins = slice_fr_interval.bins_by_step(SamplingPoint(sm, stamp=sm.stamp_from_t(timebin * qu.ms)),
                                                                SamplingPoint(sm, t=timestep * qu.ms))
            timebin_fr_dict[event_entry.Index] = [session_spiketrain.slice_by_index(index_interval).firing_rate()
                                                  for index_interval in slice_fr_time_bins]


    md.np_saver(timebin_fr_dict, target_filename)


main()
