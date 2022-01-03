import sys
import re
from ast import literal_eval
from operator import itemgetter
from itertools import compress
from scipy.io import loadmat
from metadata import *

'''Prereqs
sessions.csv: manually update with days to be considered
anatomy.csv: obtained from parseAnatomy
sorted.txt: ls -1 /jukebox/buschman/Projects/GatingInWorkingMemory/Data/*/18*/CellSorting/*.mat | xargs -n 1 basename
trialConditions.csv: obtained from parseTrialConditions(sessions)
sorting.csv: obtained from parseUnits(sessions,cellSorting)
need units entry to have a corresponding sorted .mat file
some things are saved based on their index (0-based: StageIndex, StampIndex)
and some others based on enumeration (1-based: TrialNum, ChanNum)

Argument list:
[fsroot]
'''


def main():

    # TODO: obtain event.StageDuration from bhv files Stage structs

    args = sys.argv
    
    # INIT
    md = MetaData()

    # IMPORT
    # import sessions
    sessions = pd.read_csv(md.preproc_src_path('sessions.csv'))
    # typecasting
    columns = ['Subject', 'Date', 'Sorter', 'Paused', 'Sleep', 'EncodeBug', 'BitcodeBug', 'NEVCut', 'NS4Cut', 'Sort', 'Valid']
    types = ['category', 'str', 'category', 'bool', 'bool', 'bool', 'bool', 'bool', 'bool', 'bool', 'bool']
    sessions = sessions.astype(dict(zip(columns, types)))
    sessions['DateStamp'] = pd.to_datetime(sessions['Date'], format='%y%m%d')
    sessions.Sleep = False
    # set index
    sessions.set_index(md.preproc_imports['sessions']['index'], inplace=True, drop=False)

    # import trial conditions
    trials = pd.read_csv(md.preproc_src_path('trialConditions.csv'))
    trials['Next'] = trials['Next'].fillna(-1).astype('int64')
    trials['Prev'] = trials['Prev'].fillna(-1).astype('int64')
    trials.rename({'Index': 'TrialIndex', 'Prev': 'TrialPrevIndex', 'Next': 'TrialNextIndex'}, inplace=True)
    trials['Response'] = trials['Response'].fillna(-1)
    trials['RuleCueCategory'] = 'C' + trials['RuleGroup'].astype(str) + trials['RuleCue'].astype(str)
    trials['RuleStimCategory'] = 'S' + trials['RuleGroup'].astype(str) + trials['RuleStim'].astype(str)
    # typecasting
    columns = ['Session', 'Cue', 'Stim', 'RuleGroup', 'RuleCue', 'RuleStim', 'StopCondition', 'Response', 'TargetStim', 'RuleCueCategory', 'RuleStimCategory']
    types = ['category', 'category', 'category', 'category', 'category', 'category', 'category', 'int64', 'category', 'category', 'category']
    trials = trials.astype(dict(zip(columns, types)))
    # set index
    trials.set_index(md.preproc_imports['trials']['index'], inplace=True, drop=False)

    # import and process events
    with open(md.preproc_src_path('trialCutting.txt')) as fid:
        cut_info_str = fid.read().splitlines()
    cut_info_strlist = [re.split(' ', x) for x in cut_info_str]
    cut_info_list = [[literal_eval(y) for y in x] for x in cut_info_strlist]
    cut_info_array = np.reshape(cut_info_list, [5, -1], order='F')
    cut_info_array[0] = [x[0] for x in cut_info_array[0]]
    def seq_to_frame(l, type='void'):
        if any(v is None for v in l):
            return np.tile(np.nan, len(l))
        elif type == 'stamp':
            return np.array(l) - 1
        else:
            return np.array(l)
    cutting = pd.DataFrame()
    cutting['SequenceStage'] = [seq_to_frame(l) for l in cut_info_array[1]]
    cutting['SequenceStim'] = [seq_to_frame(l) for l in cut_info_array[2]]
    cutting['SequenceEncodeIndex'] = [seq_to_frame(l) for l in cut_info_array[3]]
    cutting['SequenceBitcodeIndex'] = [seq_to_frame(l) for l in cut_info_array[4]]
    # convert to events
    events = pd.DataFrame([(trials.iloc[t].Session, trials.iloc[t].TrialNum, st,
                            cutting.iloc[t].SequenceStage[st], cutting.iloc[t].SequenceStim[st],
                            cutting.iloc[t].SequenceEncodeIndex[st], cutting.iloc[t].SequenceBitcodeIndex[st])
                           for t in range(len(trials))
                           for st in range(len(cutting.iloc[t].SequenceStage))],
                          columns=['Session', 'TrialNum', 'StageIndex', 'StageCode', 'StageStim',
                                   'StageEncodeIndex', 'StageBitcodeIndex'])
    # stagestim recode
    stage_code = list(events.StageCode)
    stage_stim = list(events.StageStim)
    stage_stim_category = [md.image_recode(st_code, im)
                           for st_code, im
                           in zip(stage_code, stage_stim)]
    events['StageStimCategory'] = stage_stim_category
    # delay of stim
    prev_stage_stim_category = [-1] + stage_stim_category[:-1]
    events['StageStimDelayCategory'] = [prev_stim if (st_code == md.specs['DelayOnset'])
                                        else 'NODELAY' if (st_code in [md.specs['FixOnset'], md.specs['CueOnset'], md.specs['StimOnset']])
                                        else None
                                        for prev_stim, st_code
                                        in zip(prev_stage_stim_category, stage_code)]
    # fix, cue, cuedelay, stim, stimdelay encoding
    prev_stage_code = [-1] + stage_code[:-1]
    events['StageCategory'] = [md.stage_to_stage_category(stage, prev_stage)
                               for stage, prev_stage
                               in zip(stage_code, prev_stage_code)]

    # estimate stage bitcode-based index
    bitcode_offset = round(np.nanmean(np.array(events.StageBitcodeIndex) - np.array(events.StageEncodeIndex)))
    encode_index = list(events.StageEncodeIndex)
    bitcode_index = list(events.StageBitcodeIndex)
    inferred_bitcode_index = list(np.array(encode_index)+bitcode_offset)
    missing_bitcode_index = [np.isnan(bitcode) and not(np.isnan(encode))
                             for encode, bitcode
                             in zip(encode_index, bitcode_index)]
    events['StageTimeIndex'] = [inferred_bitcode if missing_bitcode else bitcode
                                for inferred_bitcode, bitcode, missing_bitcode
                                in zip(inferred_bitcode_index, bitcode_index, missing_bitcode_index)]
    # stage duration in timestamps
    continues_from_prev_event = [not(stage == 2) for stage in stage_code]
    continues_on_next_event = continues_from_prev_event[1:]+[False]
    next_time_index = list(events.StageTimeIndex)[1:]+[np.nan]
    events['StageDuration'] = [next_time-time if continues else np.nan
                               for time, next_time, continues
                               in zip(list(events.StageTimeIndex), next_time_index, continues_on_next_event)]

    # typecasting
    columns = ['Session', 'StageCode', 'StageStim', 'StageStimCategory', 'StageStimDelayCategory', 'StageCategory']
    types = ['category', 'category', 'category', 'category', 'category', 'category']
    events = events.astype(dict(zip(columns, types)))
    # set index
    events.set_index(md.preproc_imports['events']['index'], inplace=True, drop=False)

    # import units
    # find sorted channels
    with open(md.preproc_src_path('sorted.txt')) as fid:
        sorted_filenames = fid.read().splitlines()
    # parse filename
    sorted_parts = list(map(lambda part: itemgetter(1, 2, 3, 5)(re.findall('[a-z]+|[0-9]+', part)), sorted_filenames))
    sorted_parts = pd.DataFrame(sorted_parts, columns=['Subject', 'Date', 'SessID', 'ChanNum'])
    sorted_parts.Subject = sorted_parts.Subject.apply(lambda name: name.capitalize())
    # construct Session string
    sorted_parts['Session'] = list(map(lambda subject, date, sessionid: "{0}_{1}_{2}".format(subject, date, sessionid),
                                       sorted_parts.Subject, sorted_parts.Date, sorted_parts.SessID))
    # remove unnecessary columns
    sorted_parts.drop(['Subject', 'Date', 'SessID'], axis=1, inplace=True)
    # typecasting
    columns = ['Session', 'ChanNum']
    types = ['category', 'int']
    sorted_parts = sorted_parts.astype(dict(zip(columns, types)))
    # import anatomy
    anatomy = pd.read_csv(md.preproc_src_path('anatomy.csv'))
    # typecasting
    columns = ['Session', 'Area', 'ChanType', 'ChanGrp']
    types = ['category', 'category', 'category', 'category']
    anatomy = anatomy.astype(dict(zip(columns, types)))
    # import sorting
    sorting = pd.read_csv(md.preproc_src_path('sorting.csv'))
    # typecasting
    sorting['RatingCode'] = sorting['RatingCode'].fillna(-1).astype('int64')
    columns = ['Session']
    types = ['category']
    sorting = sorting.astype(dict(zip(columns, types)))
    # merge dataframes for final units dataframe and set index
    units = pd.merge(pd.merge(anatomy, sorted_parts), sorting).set_index(md.preproc_imports['units']['index'], drop=False)

    physiology = pd.DataFrame(index=units.index.copy()).reset_index()

    channel_ind_cols = md.preproc_imports['multiunits']['index']
    multiunits = pd.merge(units.reset_index(drop=True).groupby(channel_ind_cols)['Session', 'ChanNum', 'Area', 'ChanType', 'ChanGrp'].head(1),
                          units.reset_index(drop=True).groupby(channel_ind_cols)['NumSpikes'].sum().reset_index()).set_index(channel_ind_cols, drop=False)



    #import activity
    unit_stats = []
    for ui in range(len(units)):
        session = units.iloc[ui].Session
        subject = sessions.loc[session].Subject
        date = sessions.loc[session].Date
        sessid = sessions.loc[session].SessID
        channum = units.iloc[ui].ChanNum
        unitnum = units.iloc[ui].UnitNum
        sortmat = loadmat(md.preproc_dat_spike_path(subject, str(date),
                                                    'gatingwm_{0}_{1}_{2:03d}_chan{3:03d}.mat'.
                                                    format(subject.lower(), date, sessid, channum)))

        def flatten(listlist):
            return [item for innerlist in listlist for item in innerlist]

        unit_inds = np.array(flatten(sortmat['units'] == unitnum))
        ts = np.array(flatten(sortmat['ts'][unit_inds]))
        # parse timestamps to segments of active time
        split_points = list(compress(range(len(ts) - 1), np.diff(ts) / 30000 > md.active_unit_window))
        segment_split = [ts] if not split_points else np.array_split(ts, np.array(split_points) + 1)
        unit_stats.extend([(session, channum, unitnum,
                            seg[0], seg[-1], seg[-1] - seg[0], len(seg),
                            (0 if not seg[-1] - seg[0] else 30e3 * len(seg) / (seg[-1] - seg[0])))
                           for seg in segment_split])

    activity = pd.DataFrame(unit_stats, columns=md.preproc_imports['units']['index'] +
                            ['SegmentStart', 'SegmentEnd', 'SegmentDuration', 'SegmentNumSpikes', 'SegmentFiringRate'])
    # typecasting
    columns = ['Session']
    types = ['category']
    activity = activity.astype(dict(zip(columns, types)))
    activity.set_index(md.preproc_imports['activity']['index'], inplace=True, drop=False)



    # import LFP

    # SAVE
    md.db_base_saver([('sessions', sessions), ('trials', trials), ('events', events), ('units', units), ('physiology', physiology), ('multiunits', multiunits), ('activity', activity)])


main()

# fr = ((activity['SegmentDuration'] * activity['SegmentFiringRate']).groupby(level=[0, 1, 2]).sum()).\
#     divide((activity['SegmentDuration'].groupby(level=[0, 1, 2]).sum()))
