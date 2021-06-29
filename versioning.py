import pandas as pd
from rec_format import interaction_term
from itertools import product, combinations, chain


# PARSER #

def parse_var(s):
    """ Parse a key, value pair, separated by '=' That's the reverse of ShellArgs.
    On the command line (argparse) a declaration will typically look like:
    foo=hello or foo="hello world" """
    items = s.split('=')
    key = items[0].strip()  # we remove blanks around keys, as is logical
    if len(items) > 1:
        # rejoin the rest:
        value = '='.join(items[1:])
        return key, value
    else:
        return None, None


def parse_vars(items):
    """Parse a series of key-value pairs and return a dictionary """
    d = {}
    if items:
        for item in items:
            key, value = parse_var(item)
            d[key] = value
    return d


def job_scheduler(args_version, args_from_parse_func=None):

    parse_version = parse_vars(args_version)
    if 'job_id' in parse_version.keys():
        return parse_vars(args_from_parse_func(parse_version))
    else:
        return parse_version


# UTILS #

def filter_df_wrapper(df, column, wrapper, arg):

    mask = wrapper(df[column], arg)
    return {'mask': mask,
            'df': df.loc[mask]}


# VERSIONING #

def version_fr_params(version_fr):

    if version_fr == 'Concat':
        t_start = -50
        t_end = 1000
        timebin = 50
        timestep = 25
        event_mask = {'column': 'StageCategory', 'wrapper': pd.Series.isin, 'arg': ['CueOnset', 'StimOnset']}
    elif version_fr == 'ConcatFactor':
        t_start = -300
        t_end = 1300
        timebin = 150
        timestep = 25
        event_mask = {'column': 'StageCategory', 'wrapper': pd.Series.isin, 'arg': ['CueOnset', 'StimOnset']}
    elif version_fr == 'ConcatFactor2':
        t_start = -600
        t_end = 1600
        timebin = 150
        timestep = 25
        event_mask = {'column': 'StageCategory', 'wrapper': pd.Series.isin, 'arg': ['CueOnset', 'StimOnset']}
    elif version_fr == 'ConcatExtended':
        t_start = -50
        t_end = 2000
        timebin = 50
        timestep = 25
        event_mask = {'column': 'StageCategory', 'wrapper': pd.Series.isin, 'arg': ['CueOnset', 'StimOnset']}
    elif version_fr == 'ConcatFactorExtended':
        t_start = -400
        t_end = 2400
        timebin = 150
        timestep = 25
        event_mask = {'column': 'StageCategory', 'wrapper': pd.Series.isin, 'arg': ['CueOnset', 'StimOnset']}
    elif version_fr == 'Sample':
        t_start = -50
        t_end = 400
        timebin = 50
        timestep = 25
        event_mask = {'column': 'StageCategory', 'wrapper': pd.Series.isin, 'arg': ['CueOnset', 'StimOnset']}
    elif version_fr == 'Delay':
        t_start = -50
        t_end = 600
        timebin = 50
        timestep = 25
        event_mask = {'column': 'StageCategory', 'wrapper': pd.Series.isin, 'arg': ['CueDelay', 'StimDelay']}
    elif version_fr == 'WindowFix':
        t_start = -500
        t_end = 0
        timebin = t_end - t_start
        timestep = timebin
        event_mask = {'column': 'StageCategory', 'wrapper': pd.Series.isin, 'arg': ['CueOnset']}
    elif version_fr == 'WindowSample':
        t_start = 0
        t_end = 385
        timebin = t_end - t_start
        timestep = timebin
        event_mask = {'column': 'StageCategory', 'wrapper': pd.Series.isin, 'arg': ['CueOnset', 'StimOnset']}
    elif version_fr == 'WindowDelay':
        t_start = 0
        t_end = 565
        timebin = t_end - t_start
        timestep = timebin
        event_mask = {'column': 'StageCategory', 'wrapper': pd.Series.isin, 'arg': ['CueDelay', 'StimDelay']}
    elif version_fr == 'WindowSampleShift':
        t_start = 50
        t_end = 435
        timebin = t_end - t_start
        timestep = timebin
        event_mask = {'column': 'StageCategory', 'wrapper': pd.Series.isin, 'arg': ['CueOnset', 'StimOnset']}
    elif version_fr == 'WindowDelayShift':
        t_start = 50
        t_end = 615
        timebin = t_end - t_start
        timestep = timebin
        event_mask = {'column': 'StageCategory', 'wrapper': pd.Series.isin, 'arg': ['CueDelay', 'StimDelay']}
    elif version_fr == 'WindowDelayLate':
        t_start = 200
        t_end = 565
        timebin = t_end - t_start
        timestep = timebin
        event_mask = {'column': 'StageCategory', 'wrapper': pd.Series.isin, 'arg': ['CueDelay', 'StimDelay']}
    elif version_fr == 'WindowGatingClassify':
        t_start = 200
        t_end = 500
        timebin = t_end - t_start
        timestep = timebin
        event_mask = {'column': 'StageCategory', 'wrapper': pd.Series.isin, 'arg': ['CueOnset', 'StimOnset']}
    elif version_fr == 'WindowMemoryClassify':
        t_start = 600
        t_end = 900
        timebin = t_end - t_start
        timestep = timebin
        event_mask = {'column': 'StageCategory', 'wrapper': pd.Series.isin, 'arg': ['CueOnset', 'StimOnset']}

    return {'t_start': t_start,
            't_end': t_end,
            'timebin': timebin,
            'timestep': timestep,
            'event_mask': event_mask}


def version_units_area_list(version_units):

    area_list = version_units.split('_')
    area_list_str = ''.join(sorted(area_list))

    return {'area_list': area_list,
            'area_list_str': area_list_str}


def version_units_subject_list(version_subjects):

    subject_list = version_subjects.split('_')
    subject_list_str = ''.join(sorted(subject_list))

    return {'subject_list': subject_list,
            'subject_list_str': subject_list_str}


# ANOVA #

def anova_version_aov_params(version_aov, version_fr):

    if version_aov == 'GeneralizedGating':
        selection_dict = {'column': 'Catch',
                          'list': ['All']}
        x_a = 'GatingCondSpecialized'
        x_b = 'StageStimSpecialized'
        x_ab = interaction_term(x_a, x_b)
        levels_dict = {x_a: ['Gating', 'PreDist'],
                       x_b: ['S11', 'S12', 'S21', 'S22']}
        if 'Delay' in version_fr:
            levels_dict = {x_a: ['GatingDelay', 'PreDistDelay'],
                           x_b: ['S11D', 'S12D', 'S21D', 'S22D']}
        event_cnj_mask = [{'column': 'StageStimSpecialized',
                           'wrapper': pd.Series.ne,
                           'arg': 'S00'}]
        if 'Delay' in version_fr:
            event_cnj_mask = [{'column': 'StageStimSpecialized',
                               'wrapper': pd.Series.ne,
                               'arg': 'S00D'}]
        group_column_list = ['RuleGroup']
        x_factors = [x_a, x_b, x_ab]
    elif version_aov == 'GeneralizedOutputGating':
        selection_dict = {'column': 'Catch',
                          'list': ['All']}
        x_a = 'GatingCondSpecialized'
        x_b = 'StageStimSpecialized'
        x_ab = interaction_term(x_a, x_b)
        levels_dict = {x_a: ['Target', 'PostDist'],
                       x_b: ['S11', 'S12', 'S21', 'S22']}
        if 'Delay' in version_fr:
            levels_dict = {x_a: ['TargetDelay', 'PostDistDelay'],
                           x_b: ['S11D', 'S12D', 'S21D', 'S22D']}
        event_cnj_mask = [{'column': 'StageStimSpecialized',
                           'wrapper': pd.Series.ne,
                           'arg': 'S00'}]
        if 'Delay' in version_fr:
            event_cnj_mask = [{'column': 'StageStimSpecialized',
                               'wrapper': pd.Series.ne,
                               'arg': 'S00D'}]
        group_column_list = ['RuleGroup']
        x_factors = [x_a, x_b, x_ab]
    elif version_aov == 'PresentedStimulus':
        selection_dict = {'column': 'GatingCondSpecialized',
                          'list': ['PreDist', 'Gating', 'PostDist', 'Target']}
        x_a = 'StageStimSpecialized'
        levels_dict = {x_a: ['S11', 'S12', 'S21', 'S22']}
        event_cnj_mask = [{'column': 'StageStimSpecialized',
                           'wrapper': pd.Series.ne,
                           'arg': 'S00'}]
        group_column_list = ['RuleGroup']
        x_factors = [x_a]
    elif version_aov == 'PresentedStimulusExtended':
        selection_dict = {'column': 'GatingCondSpecialized',
                          'list': ['PreDist', 'Gating', 'PostDist']}
        x_a = 'StageStimExtendedNoTarget'
        levels_dict = {x_a: ['S11', 'S12', 'S21', 'S22']}
        event_cnj_mask = [{'column': 'StageStimExtendedNoTarget',
                           'wrapper': pd.Series.ne,
                           'arg': 'S00'}]
        group_column_list = ['RuleGroup']
        x_factors = [x_a]
    elif version_aov == 'NextStimulusExtended':
        selection_dict = {'column': 'GatingCondSpecialized',
                          'list': ['PreDist', 'Gating', 'PostDist']}
        x_a = 'NextStageStimExtendedNoTarget'
        levels_dict = {x_a: ['S11', 'S12', 'S21', 'S22']}
        event_cnj_mask = [{'column': 'NextStageStimExtendedNoTarget',
                           'wrapper': pd.Series.ne,
                           'arg': 'S00'}]
        group_column_list = ['RuleGroup']
        x_factors = [x_a]
    elif version_aov == 'GatedStimulus':
        selection_dict = {'column': 'GatingCondSpecialized',
                          'list': ['Cue', 'PreDist', 'Gating', 'PostDist', 'Target']}
        x_a = 'RuleStimCategory'
        levels_dict = {x_a: ['S11', 'S12', 'S21', 'S22']}
        event_cnj_mask = [{'column': 'StageStimSpecialized',
                           'wrapper': pd.Series.ne,
                           'arg': 'S00'}]
        group_column_list = ['RuleGroup']
        x_factors = [x_a]
    elif version_aov == 'GatedStimulusPostDistractor':
        selection_dict = {'column': 'GatedStimulusSerialPosition',
                          'list': ['Gating', 'PostDist1', 'PostDist2', 'Target']}
        if 'Delay' in version_fr:
            selection_dict = {'column': 'GatedStimulusSerialPosition',
                              'list': ['Gating', 'PostDist1', 'PostDist2']}
        x_a = 'RuleStimCategory'
        levels_dict = {x_a: ['S11', 'S12', 'S21', 'S22']}
        event_cnj_mask = [{'column': 'StageStimSpecialized',
                           'wrapper': pd.Series.ne,
                           'arg': 'S00'}]
        if 'Delay' in version_fr:
            event_cnj_mask = [{'column': 'StageStimSpecialized',
                               'wrapper': pd.Series.ne,
                               'arg': 'S00D'}]
        group_column_list = ['RuleGroup']
        x_factors = [x_a]
    elif version_aov == 'GatedStimulusPostDistractorCentered':
        selection_dict = {'column': 'GatedStimulusSerialPositionCentered',
                          'list': ['Gating-1', 'Gating', 'Gating+1', 'Gating+2']}
        x_a = 'RuleStimCategory'
        levels_dict = {x_a: ['S11', 'S12', 'S21', 'S22']}
        event_cnj_mask = [{'column': 'StageStimSpecialized',
                           'wrapper': pd.Series.ne,
                           'arg': 'S00'}]
        if 'Delay' in version_fr:
            event_cnj_mask = [{'column': 'StageStimSpecialized',
                               'wrapper': pd.Series.ne,
                               'arg': 'S00D'}]
        group_column_list = ['RuleGroup']
        x_factors = [x_a]
    elif version_aov == 'GatedGroup':
        selection_dict = {'column': 'GatingCondSpecialized',
                          'list': ['Cue', 'PreDist', 'Gating', 'PostDist', 'Target']}
        x_a = 'RuleGroup'
        levels_dict = {x_a: [1, 2]}
        event_cnj_mask = [{'column': 'StageStimSpecialized',
                           'wrapper': pd.Series.ne,
                           'arg': 'S00'}]
        group_column_list = []
        x_factors = [x_a]
    elif version_aov == 'NullGatedStimulus':
        selection_dict = {'column': 'GatingCondSpecialized',
                          'list': ['PreDist', 'PostDist']}
        x_a = 'RuleStimCategory'
        levels_dict = {x_a: ['S11', 'S12', 'S21', 'S22']}
        event_cnj_mask = [{'column': 'StageStimSpecialized',
                           'wrapper': pd.Series.eq,
                           'arg': 'S00'}]
        group_column_list = ['RuleGroup']
        x_factors = [x_a]
    elif version_aov == 'NullGatedGroup':
        selection_dict = {'column': 'GatingCondSpecialized',
                          'list': ['PreDist', 'PostDist']}
        x_a = 'RuleGroup'
        levels_dict = {x_a: [1, 2]}
        event_cnj_mask = [{'column': 'StageStimSpecialized',
                           'wrapper': pd.Series.eq,
                           'arg': 'S00'}]
        group_column_list = []
        x_factors = [x_a]
    elif version_aov == 'PreviousStimulus':
        selection_dict = {'column': 'GatingCond_From_To_GatingCondSpecialized',
                          'list': ['PreDist_From_To_PreDist', 'PreDist_From_To_Gating', 'Gating_From_To_PostDist',
                                   'PostDist_From_To_PostDist', 'PostDist_From_To_Target']}
        x_a = 'PrevStageStimSpecialized'
        levels_dict = {x_a: ['S11', 'S12', 'S21', 'S22']}
        event_cnj_mask = [{'column': 'PrevStageStimSpecialized',
                           'wrapper': pd.Series.ne,
                           'arg': 'S00'}]
        group_column_list = ['RuleGroup']
        x_factors = [x_a]
    elif version_aov == 'PreviousDistractor':
        selection_dict = {'column': 'Catch',
                          'list': ['All']}
        x_a = 'PrevDistractorSpecialized'
        levels_dict = {x_a: ['S11', 'S12', 'S21', 'S22']}
        event_cnj_mask = [{'column': 'PrevDistractorSpecialized',
                           'wrapper': pd.Series.ne,
                           'arg': 'S00'}]
        group_column_list = ['RuleGroup']
        x_factors = [x_a]

    return {'selection_dict': selection_dict,
            'levels_dict': levels_dict,
            'event_cnj_mask': event_cnj_mask,
            'group_column_list': group_column_list,
            'x_factors': x_factors}


# FACTOR / DPCA #

def factor_generate_conditions(version_factor):

    if version_factor == 'StimulusGating':

        condition_columns = ['StageStimSpecialized', 'GatingCondSpecialized']
        condition_list = list(product(['S11', 'S12', 'S21', 'S22'], ['PreDist', 'Gating', 'PostDist', 'Target']))
        balance_columns = []

    elif version_factor == 'StimulusGatingBool':

        condition_columns = ['StageStimSpecialized', 'GatingBoolSpecialized']
        condition_list = list(product(['S11', 'S12', 'S21', 'S22'], ['Gating', 'Dist']))
        balance_columns = []

    elif version_factor == 'StimulusGatingPreBool':

        condition_columns = ['StageStimSpecialized', 'GatingCondSpecialized']
        condition_list = list(product(['S11', 'S12', 'S21', 'S22'], ['Gating', 'PreDist']))
        balance_columns = []

    elif version_factor == 'RuleStimGating':

        condition_columns = ['RuleStimCategory', 'GatingCondSpecialized']
        condition_list = list(product(['S11', 'S12', 'S21', 'S22'], ['Cue', 'PreDist', 'Gating', 'PostDist', 'Target']))
        balance_columns = []

    elif version_factor == 'RuleStimGatingBool':

        condition_columns = ['RuleStimCategory', 'GatingBoolSpecialized']
        condition_list = list(product(['S11', 'S12', 'S21', 'S22'], ['Gating', 'Dist']))
        balance_columns = []

    elif version_factor == 'PostStimulusRuleStim':

        condition_columns = ['PostStageStimSpecialized', 'PostRuleStimCategory']
        condition_list = list(product(['S11', 'S12', 'S21', 'S22'], ['S11', 'S12', 'S21', 'S22']))
        balance_columns = []

    elif version_factor == 'GatPostStimulusRuleStim':

        condition_columns = ['GatPostStageStimSpecialized', 'GatPostRuleStimCategory']
        condition_list = list(product(['S11', 'S12', 'S21', 'S22'], ['S11', 'S12', 'S21', 'S22']))
        balance_columns = []

    elif version_factor == 'RuleStimGatingNull':

        condition_columns = ['RuleStimCategory', 'GatingNullCondSpecialized']
        condition_list = list(product(['S11', 'S12', 'S21', 'S22'], ['PreDist', 'PostDist']))
        balance_columns = []

    elif version_factor == 'GatedStimulus':

        condition_columns = ['GatedStageStimSpecialized']
        condition_list = ['S11', 'S12', 'S21', 'S22']
        balance_columns = []

    elif version_factor == 'GatedStimulusPostDistMemory':

        condition_columns = ['GatedStimulusPostDistMemory', 'SM_SensoryAbstractGroup']
        condition_list = list(product(['S11', 'S12', 'S21', 'S22'], ['AbstractGroup1', 'AbstractGroup2']))
        balance_columns = ['SM_SensoryAbstractGroup']

    elif version_factor == 'GatingPreBool':

        condition_columns = ['GatingCondSpecialized', 'StageStimSpecialized']
        condition_list = list(product(['Gating', 'PreDist'], ['S11', 'S12', 'S21', 'S22']))
        balance_columns = ['StageStimSpecialized']


    return {'condition_columns': condition_columns,
            'condition_list': condition_list,
            'balance_columns': balance_columns}


def factor_version_filter_params(version_filter):

    filter_params_list = version_filter.split('_')
    filter_params_list_str = ''.join(sorted(filter_params_list))

    balance = False
    average = False

    if 'Balance' in filter_params_list:
        balance = True
    if 'Average' in filter_params_list:
        average = True

    return {'balance': balance,
            'average': average,
            'filter_params_list_str': filter_params_list_str}


def factor_dpca_labels_mapping(version_factor):

    if version_factor == 'StimulusGating':
        return 'sgt'
    elif version_factor == 'RuleStimGating':
        return 'mgt'
    elif version_factor == 'PostStimulusRuleStim':
        return 'smt'
    elif version_factor == 'GatPostStimulusRuleStim':
        return 'smt'
    elif version_factor == 'StimulusGatingBool':
        return 'sgt'
    elif version_factor == 'StimulusGatingPreBool':
        return 'sgt'
    elif version_factor == 'RuleStimGatingBool':
        return 'mgt'
    elif version_factor == 'RuleStimGatingNull':
        return 'mgt'
    elif version_factor == 'GatedStimulus':
        return 'st'
    elif version_factor == 'GatedStimulusPostDistMemory':
        return 'mt'
    elif version_factor == 'GatingPreBool':
        return 'gt'


def factor_dpca_join_mapping(labels):

    cond_labels = labels.replace('t', '')
    num_margins = len(cond_labels)
    combs = list(chain.from_iterable(combinations(list(range(num_margins)), r + 1) for r in range(num_margins + 1)))
    def get_letters_from_inds(inds): return ''.join([cond_labels[ind] for ind in inds]) if bool(inds) else 't'
    def letters_list(letters): return [letters, letters + 't'] if letters != 't' else ['t']
    return {get_letters_from_inds(inds): letters_list(get_letters_from_inds(inds)) for inds in combs}


# CLASSIFICATION #

def classification_version_class(version_class):
    """ Defines variable to classify """

    if version_class == 'GatingPreBool':

        condition_columns = ['GatingCondSpecialized']
        condition_list = [['PreDist', 'Gating']]

    elif version_class == 'Stimulus':

        condition_columns = ['StageStimSpecialized']
        condition_list = [['S11', 'S12', 'S21', 'S22']]

    elif version_class == 'GatedStimulus':

        condition_columns = ['RuleStimCategory']
        condition_list = [['S11', 'S12', 'S21', 'S22']]

    return {'condition_columns': condition_columns,
            'condition_list': condition_list}


def classification_version_balance(version_balance):

    if version_balance == 'Stimulus':

        condition_columns = ['StageStimSpecialized']
        condition_list = [['S11', 'S12', 'S21', 'S22']]
        split_ignore_columns = []

    elif version_balance == 'StageGatingPrePost':

        condition_columns = ['GatingCondSpecialized']
        condition_list = [['PreDist', 'Gating', 'PostDist']]
        split_ignore_columns = []

    elif version_balance == 'StageGatingCentered':

        condition_columns = ['GatedStimulusSerialPositionCentered']
        condition_list = [['Gating-1', 'Gating', 'Gating+1']]
        split_ignore_columns = []

    elif version_balance == 'StageGatingPrePostSensory':

        condition_columns = ['GatingCondSpecialized', 'SM_SensoryAbstractGroup']
        condition_list = [['PreDist', 'Gating', 'PostDist'], ['AbstractGroup1', 'AbstractGroup2']]
        split_ignore_columns = ['SM_SensoryAbstractGroup']

    elif version_balance == 'StageGatingCenteredSensory':

        condition_columns = ['GatedStimulusSerialPositionCentered', 'SM_SensoryAbstractGroup']
        condition_list = [['Gating-1', 'Gating', 'Gating+1'], ['AbstractGroup1', 'AbstractGroup2']]
        split_ignore_columns = ['SM_SensoryAbstractGroup']

    elif version_balance == 'StageGatingPrePostMemory':

        condition_columns = ['GatingCondSpecialized', 'SM_MemoryAbstractGroup']
        condition_list = [['PreDist', 'Gating', 'PostDist'], ['AbstractGroup1', 'AbstractGroup2']]
        split_ignore_columns = ['SM_MemoryAbstractGroup']

    elif version_balance == 'StageGatingCenteredMemory':

        condition_columns = ['GatedStimulusSerialPositionCentered', 'SM_MemoryAbstractGroup']
        condition_list = [['Gating-1', 'Gating', 'Gating+1'], ['AbstractGroup1', 'AbstractGroup2']]
        split_ignore_columns = ['SM_MemoryAbstractGroup']

    elif version_balance == 'StageGatingCenteredGatingOnly':

        condition_columns = ['GatedStimulusSerialPositionCentered']
        condition_list = [['Gating']]
        split_ignore_columns = []

    elif version_balance == 'StageGatingCenteredPostDist1Only':

        condition_columns = ['GatedStimulusSerialPositionCentered']
        condition_list = [['Gating+1']]
        split_ignore_columns = []

    elif version_balance == 'StageGatingCenteredMemoryGatingOnly':

        condition_columns = ['GatedStimulusSerialPositionCentered', 'SM_MemoryAbstractGroup']
        condition_list = [['Gating'], ['AbstractGroup1', 'AbstractGroup2']]
        split_ignore_columns = ['SM_MemoryAbstractGroup']

    elif version_balance == 'StageGatingCenteredMemoryPostDist1Only':

        condition_columns = ['GatedStimulusSerialPositionCentered', 'SM_MemoryAbstractGroup']
        condition_list = [['Gating+1'], ['AbstractGroup1', 'AbstractGroup2']]
        split_ignore_columns = ['SM_MemoryAbstractGroup']

    return {'condition_columns': condition_columns,
            'condition_list': condition_list,
            'split_ignore_columns': split_ignore_columns}
