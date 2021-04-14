import numpy as np
import pandas as pd
import pickle
from re import findall
from os import makedirs
from environment import *
from scipy import signal


class MetaData:

    # TODO: make pathing to different superstructs more robust and flexible

    env = Environment()

    # file system
    preproc_src_superstruct = path.join(env.preproc_src_env, 'Analysis/PreProcessing')
    preproc_dat_superstruct = env.preproc_dat_env
    preproc_dest_superstruct = path.join(env.preproc_dest_env, 'input')
    proc_dest_superstruct = path.join(env.proc_dest_env, 'output')
    preproc_imports = {'sessions': {'filename': 'sessions.pkl', 'index': ['Session']},
                       'trials': {'filename': 'trials.pkl', 'index': ['Session', 'TrialNum']},
                       'events': {'filename': 'events.pkl', 'index': ['Session', 'TrialNum', 'StageIndex']},
                       'units': {'filename': 'units.pkl', 'index': ['Session', 'ChanNum', 'UnitNum']},
                       'physiology': {'filename': 'physiology.pkl', 'index': ['Session', 'ChanNum', 'UnitNum']},
                       'multiunits': {'filename': 'multiunits.pkl', 'index': ['Session', 'ChanNum']},
                       'activity': {'filename': 'activity.pkl', 'index': ['Session', 'ChanNum', 'UnitNum']},
                       'conditions': {'filename': 'conditions.pkl', 'index': []}}
    proc_imports = {'behavioral_units': {'index': ['Session', 'ChanNum', 'UnitNum', 'TrialNum', 'StageIndex']},
                    'behavioral_multiunits': {'index': ['Session', 'ChanNum', 'TrialNum', 'StageIndex']}}

    # conditioning
    specs = {'FixOnset': 2, 'CueOnset': 11, 'DelayOnset': 10, 'StimOnset': 8}
    stopcondition_to_category = {1: 'CORRECT', -2: 'IDLE', -3: 'FIX_BREAK', -4: 'NO_BAR_HOLD', -6: 'NO_RESPONSE', -7: 'EARLY_RESPONSE'}

    # processing opts
    active_unit_window = 60*60 # measured in s
    base_rate = 30e3  # measured in Hz
    filt_deltat = 0.085
    filt_sdscale = 0.2
    filt_boxdeltat = 0.075
    filt_win_gauss = signal.windows.gaussian(filt_deltat*base_rate, filt_sdscale*filt_deltat*base_rate)
    filt_win_box = signal.windows.boxcar(int(filt_boxdeltat*base_rate))

    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 150)

    # CueOnset/StimOnset mean: 385ms
    # CueOnset/StimOnset std: 0.3ms
    # CueDelay/StimDelay mean: 567ms
    # CueDelay/StimDelay std: 0.4ms


    def preproc_src_path(self, fname):
        return path.join(self.preproc_src_superstruct, fname)

    def preproc_dat_spike_path(self, subject, date, fname):
        return path.join(self.preproc_dat_superstruct, subject, date, 'CellSorting', fname)

    def preproc_dest_path(self, fname):
        return path.join(self.preproc_dest_superstruct, fname)

    def spike_dest_path(self, datatype, session, channum, unitnum):
        return path.join(self.preproc_dest_path(path.join('spike_{0:s}'.format(datatype),
                                                          'spike_{0:s}_{1:s}_chan{2:03d}_unit{3:03d}.pkl'.
                                                          format(datatype, session, channum, unitnum))))

    def proc_dest_path(self, analysis, fname):
        dest_path = path.join(self.proc_dest_superstruct, analysis)
        if not path.exists(dest_path):
            makedirs(dest_path)
        return path.join(dest_path, fname)

    def df_saver(self, table, path):
        table.reset_index(drop=True).to_pickle(path)

    def df_loader(self, path, index):
        df = pd.read_pickle(path)
        try:
            df.set_index(index, drop=False, inplace=True)
            return df
        except (KeyError, ValueError) as e:
            return df

    def db_base_loader(self, keys_list=None):
        db = {}
        keys = self.preproc_imports.keys() if keys_list is None else keys_list
        for k in keys:
            db[k] = self.df_loader(self.preproc_dest_path(self.preproc_imports[k]['filename']), self.preproc_imports[k]['index'])
        return db

    def db_base_saver(self, df_list):
        for k, df in df_list:
            self.df_saver(df, self.preproc_dest_path(self.preproc_imports[k]['filename']))

    def np_saver(self, variable, path):
        with open(path, 'wb') as f:
            pickle.Pickler(f).dump(variable)

    def np_loader(self, path):
        return pd.read_pickle(path)
        # with open(path, 'rb') as f:
        #     return pickle.Unpickler(f).load()
    
    def image_to_group_image_pair(self, im):
        group_ind, im_ind = np.unravel_index(int(im)-1, (2, 2), order='F')
        return (group_ind+1, im_ind+1)
    
    def image_recode(self, st_code, im):
        if st_code in [self.specs['FixOnset'], self.specs['DelayOnset']]:
            return 'FIX'
        elif st_code == self.specs['CueOnset']:
            group, cue = self.image_to_group_image_pair(im)
            return 'C{0:d}{1:d}'.format(group, cue)
        elif st_code == self.specs['StimOnset']:
            if im == 0:
                return 'S00'
            else:
                group, stim = self.image_to_group_image_pair(im)
                return 'S{0:d}{1:d}'.format(group, stim)
        else:
            return None

    def image_decode(self, st_cat):
        if type(st_cat) == str and len(findall('[C|S]\d\d', st_cat)) > 0:
            return int(st_cat[1]), int(st_cat[2])
        else:
            return None

    def stage_to_stage_category(self, st_code, prev_st_code):
        if st_code == self.specs['FixOnset']:
            return 'FixOnset'
        elif st_code == self.specs['CueOnset']:
            return 'CueOnset'
        elif st_code == self.specs['StimOnset']:
            return 'StimOnset'
        elif st_code == self.specs['DelayOnset']:
            if prev_st_code == self.specs['CueOnset']:
                return 'CueDelay'
            elif prev_st_code == self.specs['StimOnset']:
                return 'StimDelay'
        else:
            return None

    def stage_to_gating_cond(self, st_ind, numpre, numpost):
        if 1 <= st_ind <= 2:
            return 'Cue'
        elif 3 <= st_ind <= 2*(numpre+1):
            return 'PreDist'
        elif 2*(numpre+1)+1 <= st_ind <= 2*(numpre+2):
            return 'Gating'
        elif 2*(numpre+2)+1 <= st_ind <= 2*(numpre+numpost+2):
            return 'PostDist'
        elif 2*(numpre+2+numpost)+1 <= st_ind <= 2*(numpre+numpost+3):
            return 'Target'

    def df_GatingCondStageStimCategory(self, row):

        if pd.isna(row['GatingCondExtended']):
            return np.nan
        if row['GatingCondExtended'] in ['Cue', 'Gating', 'Target']:
            return row['GatingCondExtended']
        elif row['StageStimExtended'] == 'S00':
            return row['GatingCondExtended'] + 'Null'
        else:
            DistGroup = 'Same' if int(row['StageStimExtended'][1]) == int(row['RuleGroup']) else 'Diff'
            return row['GatingCondExtended'] + DistGroup

    def df_BarStatus(self, row):

        valid_trial = row['StopCondition'] in [1, -3, -5, -6, -7, -8]
        attempted_trial = row['StopCondition'] in [1, -5, -6, -7, -8]
        last_event = row['StageIndex'] == row['StageNum'] - 1
        if not valid_trial:
            return 'EMPTY'
        elif not attempted_trial:
            return 'HOLD'
        elif last_event and row['StopCondition'] in [1, -5, -7, -8]:
            return 'RELEASE'
        elif (last_event and row['StopCondition'] in [-6]) or not last_event:
            return 'HOLD'

    def df_group_EnumerateStimOccurrence(self, slice):

        slice_stim_list = slice[1]['StageStimSpecialized'].value_counts().keys()
        stim_counter = dict([(stim, 0) for stim in slice_stim_list])
        occurrence_list = []
        for row in slice[1].itertuples():
            if pd.isna(row.StageStimSpecialized) or row.StageStimSpecialized[0] == 'C':
                occurrence_list.append(np.nan)
            else:
                stim_counter[row.StageStimSpecialized] += 1
                occurrence_list.append(stim_counter[row.StageStimSpecialized])

        return occurrence_list

    def df_DistractorSerialPosition(self, row):

        def stage_index_to_enumeration_suffix(index, numpre, gatingcond):
            if gatingcond == 'PreDist':
                return str(int(np.ceil(index / 2) - 1))
            elif gatingcond == 'PostDist':
                return str(int(np.ceil((index - (5 + (2 * numpre - 1))) / 2)))
            else:
                return np.nan

        return np.nan if pd.isna(row['GatingCondExtended']) else \
               stage_index_to_enumeration_suffix(row['StageIndex'], row['NumPre'], row['GatingCondExtended'])

    def df_PostStageStimSpecialized(self, row):

        sss = row['StageStimSpecialized']
        gcs = row['GatingCondSpecialized']
        return sss if gcs in ['PostDist', 'Target'] else np.nan

    def df_PostRuleStimCategory(self, row):

        rsc = row['RuleStimCategory']
        gcs = row['GatingCondSpecialized']
        return rsc if gcs in ['PostDist', 'Target'] else np.nan

    def df_GatPostStageStimSpecialized(self, row):

        sss = row['StageStimSpecialized']
        gcs = row['GatingCondSpecialized']
        return sss if gcs in ['Gating', 'PostDist'] else np.nan

    def df_GatPostRuleStimCategory(self, row):

        rsc = row['RuleStimCategory']
        gcs = row['GatingCondSpecialized']
        return rsc if gcs in ['Gating', 'PostDist'] else np.nan

    def df_GatingNullCondSpecialized(self, row):

        gcs = row['GatingCondSpecialized']
        sss = row['StageStimSpecialized']
        return gcs if sss == 'S00' else np.nan

    def beh_unit_fr_cond_col(self, y, step, segment):
        return '{0:s}_step{1:04.0f}ms_{2:04.0f}_{2:04.0f}'.format(y, step, segment[0], segment[1])

    def interval_split_to_timebins(self, interval, binsize, binstep):
        onset_list = np.linspace(interval[0], interval[1]-binsize, int(((interval[1] - interval[0] - binsize) / binstep) + 1), endpoint=True)
        return [[int(on), int(on + binsize)] for on in onset_list]