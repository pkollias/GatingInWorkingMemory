from metadata import *

'''Prereqs
MetaData class and saved sessions, trials, events, units
generates conditions table with extra columns:
...
'''


def main():

    # INIT
    md = MetaData()
    db = md.db_base_loader(['trials', 'events'])
    trials, events = db['trials'], db['events']
    trials_index = md.preproc_imports['trials']['index']
    events_index = md.preproc_imports['events']['index']
    # create conditions table
    conditions = events.reset_index(drop=True).merge(trials.reset_index(drop=True))

    # create conditions
    conditions['Catch'] = ['All' for _ in np.arange(len(conditions))]

    conditions['GatingCondExtended'] = [md.stage_to_gating_cond(st_ind, numpre, numpost)
                                for st_ind, numpre, numpost
                                in zip(list(conditions.StageIndex), list(conditions.NumPre), list(conditions.NumPost))]

    conditions['GatingCondSpecialized'] = [None if gc is None else gc + ('Delay' if sc == md.specs['DelayOnset'] else '')
                                          for gc, sc
                                          in zip(list(conditions.GatingCondExtended), list(conditions.StageCode))]

    prev_trialnum = [np.nan, np.nan] + list(events['TrialNum'])[:-2]
    prev_session = [np.nan, np.nan] + list(events['Session'])[:-2]
    prev_gatingconde = [np.nan, np.nan] + list(conditions['GatingCondExtended'])[:-2]
    prev_gatingconds = [np.nan, np.nan] + list(conditions['GatingCondSpecialized'])[:-2]
    conditions['GatingCond_From_To_GatingCondExtended'] = [np.nan if any([el == None or el == np.nan for el in [prev_gce, gce]])
                                                           else '_From_To_'.join([prev_gce, gce]) if prev_sess == sess and prev_tn == tn else np.nan
                                                           for (sess, prev_sess, tn, prev_tn, gce, prev_gce) in
                                                           zip(list(events['Session']), prev_session,
                                                               list(events['TrialNum']), prev_trialnum,
                                                               list(conditions['GatingCondExtended']), prev_gatingconde)]
    conditions['GatingCond_From_To_GatingCondSpecialized'] = [np.nan if any([el == None or el == np.nan for el in [prev_gcs, gcs]])
                                                              else '_From_To_'.join([prev_gcs, gcs]) if prev_sess == sess and prev_tn == tn else np.nan
                                                              for (sess, prev_sess, tn, prev_tn, gcs, prev_gcs) in
                                                              zip(list(events['Session']), prev_session,
                                                                  list(events['TrialNum']), prev_trialnum,
                                                                  list(conditions['GatingCondSpecialized']), prev_gatingconds)]

    conditions['StageStimExtended'] = [stim_cat if stage_cat in ['CueOnset', 'StimOnset']
                                       else stim_delay_cat if stage_cat in ['CueDelay', 'StimDelay'] else np.nan
                                       for stim_cat, stim_delay_cat, stage_cat
                                       in zip(list(conditions.StageStimCategory),
                                              list(conditions.StageStimDelayCategory),
                                              list(conditions.StageCategory))]

    conditions['StageStimSpecialized'] = [stim_cat if stage_cat in ['CueOnset', 'StimOnset']
                                       else stim_delay_cat + 'D' if stage_cat in ['CueDelay', 'StimDelay'] else np.nan
                                       for stim_cat, stim_delay_cat, stage_cat
                                       in zip(list(conditions.StageStimCategory),
                                              list(conditions.StageStimDelayCategory),
                                              list(conditions.StageCategory))]

    prev_stime = [np.nan, np.nan] + list(conditions['StageStimExtended'])[:-2]
    prev_stims = [np.nan, np.nan] + list(conditions['StageStimSpecialized'])[:-2]
    conditions['PrevStageStimExtended'] = [prev_se if prev_sess == sess and prev_tn == tn else np.nan
                                           for (sess, prev_sess, tn, prev_tn, prev_se) in
                                           zip(list(events['Session']), prev_session,
                                               list(events['TrialNum']), prev_trialnum,
                                               prev_stime)]
    conditions['PrevStageStimSpecialized'] = [prev_ss if prev_sess == sess and prev_tn == tn else np.nan
                                              for (sess, prev_sess, tn, prev_tn, prev_ss) in
                                              zip(list(events['Session']), prev_session,
                                                  list(events['TrialNum']), prev_trialnum,
                                                  prev_stims)]



    stage_stim_group = [md.image_decode(im_cat)[0] if type(im_cat) == str else np.nan
                        for im_cat in conditions.StageStimExtended]
    trial_cue_group = [md.image_to_group_image_pair(im)[0] for im in conditions.Cue]
    conditions['PostDistCategory'] = [None if gc != 'PostDist' else
                                      'Null' if ssg == 0 else 'Same' if ssg == tcg else 'Diff'
                                      for (gc, ssg, tcg) in
                                      zip(list(conditions.GatingCondExtended), stage_stim_group, trial_cue_group)]

    prev_gatingconde = [np.nan, np.nan] + list(conditions['GatingCondExtended'])[:-2]
    conditions['PrevDistractorSpecialized'] = [prev_ss if (sess == prev_sess and prev_gce in ['PreDist', 'PostDist'] and
                                                           gce in ['PreDist', 'Gating', 'PostDist'])
                                               else np.nan
                                               for (ss, sess, gce, prev_ss, prev_sess, prev_gce) in
                                               zip(list(conditions['StageStimSpecialized']),
                                                   list(events['Session']),
                                                   list(conditions['GatingCondExtended']),
                                                   prev_stims, prev_session, prev_gatingconde)]

    conditions['DistractorSerialPosition'] = conditions.apply(md.df_DistractorSerialPosition, axis=1)
    conditions['GatedStimulusSerialPosition'] = [gcs if gcs in ['Gating', 'Target']
                                                 else gcs + dsp if gcs == 'PostDist'
                                                 else np.nan
                                                 for (gcs, dsp) in
                                                 zip(list(conditions['GatingCondExtended']),
                                                     list(conditions['DistractorSerialPosition']))]
    conditions['GatingCondStageStimCategory'] = conditions.apply(md.df_GatingCondStageStimCategory, axis=1)
    conditions['BarStatus'] = conditions.apply(md.df_BarStatus, axis=1)
    occurrence_columns = events_index + ['StageStimSpecialized']
    occurrence_list = []
    for events_trial_slice in conditions[occurrence_columns].groupby(trials_index).__iter__():
        occurrence_list.extend(md.df_group_EnumerateStimOccurrence(events_trial_slice))
    conditions['StimOccurrence'] = occurrence_list

    conditions['PostStageStimSpecialized'] = conditions.apply(md.df_PostStageStimSpecialized, axis=1)
    conditions['PostRuleStimCategory'] = conditions.apply(md.df_PostRuleStimCategory, axis=1)



    # typecasting
    columns = ['Catch', 'GatingCondExtended', 'GatingCondSpecialized', 'GatingCond_From_To_GatingCondExtended', 'GatingCond_From_To_GatingCondSpecialized',
               'StageStimExtended', 'StageStimSpecialized', 'PrevStageStimExtended', 'PrevStageStimSpecialized',
               'PostDistCategory', 'RuleCueCategory', 'RuleStimCategory', 'RuleGroup', 'PrevDistractorSpecialized',
               'DistractorSerialPosition', 'GatedStimulusSerialPosition', 'GatingCondStageStimCategory', 'BarStatus', 'StimOccurrence',
               'PostStageStimSpecialized', 'PostRuleStimCategory']
    types = ['category', 'category', 'category', 'category', 'category',
             'category', 'category', 'category', 'category',
             'category', 'category', 'category', 'category', 'category',
             'category', 'category', 'category', 'category', 'category',
             'category', 'category']
    conditions = conditions.astype(dict(zip(columns, types)))
    conditions_columns = columns
    conditions = conditions[events_index + conditions_columns]

    # SAVE
    md.db_base_saver([('conditions', conditions)])


main()
