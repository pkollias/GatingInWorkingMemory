# visualize units_events

from rec_utils import *
import matplotlib.pyplot as plt

md = MetaData()
units_events = md.np_loader(md.preproc_dest_path(path.join('units_events.pkl')))
units = md.db_base_loader(['units'])['units']


sess = 'Gonzo_180120_001'
def extents(f):
  delta = f[1] - f[0]
  return [f[0] - delta/2, f[-1] + delta/2]

units_slice = units.loc[units_events[sess].columns]
units_inds = list(units_slice.loc[units_slice['Area'].eq('PFC') &
                                  units_slice['UnitNum'].ne(0) &
                                  units_slice['RatingCode'].ne(7)].index)

data = units_events[sess][units_inds].to_numpy()
x = range(data.shape[1])
y = range(data.shape[0])

plt.imshow(data, aspect='auto', interpolation='none',
           extent=extents(x) + extents(y), origin='lower')


###################

# correlation classifiers gridsearch performance

from rec_utils import *
import matplotlib.pyplot as plt

db = DataBase(['sessions', 'units', 'events'])
sessions = db.tables['sessions']

for v_ii, version_prefix_x in enumerate(['GatingPreBool_Stimulus_WindowGatingClassify_PFC',
                                         'GatingPreBool_Stimulus_WindowGatingClassify_Stri',
                                         'GatedStimulus_StageGatingCenteredSensoryGatingOnly_WindowInterferenceClassify_PFC',
                                         'GatedStimulus_StageGatingCenteredSensoryGatingOnly_WindowMemoryClassify_PFC',
                                         'GatedStimulus_StageGatingCenteredSensoryPostDist1Only_WindowInterferenceClassify_PFC',
                                         'Stimulus_StageGatingCenteredMemoryPostDist1Only_WindowInterferenceClassify_PFC',
                                         'GatedStimulus_StageGatingCenteredSensoryPostDist1Only_WindowMemoryClassify_PFC',
                                         'Stimulus_StageGatingCenteredMemoryPostDist1Only_WindowMemoryClassify_PFC']):
    version_results = []
    for sess_ind in range(42):
        print(v_ii, sess_ind)
        version_str = version_prefix_x + '_{0:d}'.format(sess_ind)
        try:
            version_str = version_prefix_x + '_{0:d}'.format(sess_ind)
            target_filename = MetaData().proc_dest_path(path.join('BehavioralUnits', 'GatingCorrelation',
                                                                  'consolidate', version_str),
                                                        'result_dict.pkl')
            result_dict = MetaData().np_loader(target_filename)
            version_results.append([res[0] for res in result_dict.values()])
        except OSError as e:
            version_results.append([])
    plt.figure()
    plt.boxplot(version_results)


###################


### correlation plots

from rec_utils import *
import matplotlib.pyplot as plt

db = DataBase(['sessions', 'units', 'events', 'conditions'])
sessions = db.tables['sessions']

try:
    max_events = MetaData().np_loader('files/max_events.pkl')
except FileNotFoundError as e:

    max_events = {}
    for v_ii, version_prefix_x in enumerate(['GatingPreBool_Stimulus_WindowGatingClassify_PFC',
                                             'GatingPreBool_Stimulus_WindowGatingClassify_Stri',
                                             'GatedStimulus_StageGatingCenteredSensoryGatingOnly_WindowMemoryClassify_PFC',
                                             'GatedStimulus_StageGatingCenteredSensoryPostDist1Only_WindowMemoryClassify_PFC',
                                             'GatedStimulus_StageGatingCenteredSensoryGatingOnly_WindowInterferenceClassify_PFC',
                                             'GatedStimulus_StageGatingCenteredSensoryPostDist1Only_WindowInterferenceClassify_PFC',
                                             'Stimulus_StageGatingCenteredMemoryPostDist1Only_WindowMemoryClassify_PFC',
                                             'Stimulus_StageGatingCenteredMemoryPostDist1Only_WindowInterferenceClassify_PFC']):
        max_events_version = []
        for sess_ind in range(42):
            print(v_ii, sess_ind)
            version_str = version_prefix_x + '_{0:d}'.format(sess_ind)
            try:
                target_filename = MetaData().proc_dest_path(path.join('BehavioralUnits', 'GatingCorrelation',
                                                                      'consolidate', version_str),
                                                            'result_dict.pkl')
                result_dict = MetaData().np_loader(target_filename)
                sorted_classifiers = sorted(result_dict.items(), key=lambda kv: 0 if kv[0][5] == 'svm' else kv[1][0])
                max_events_version.append(sorted_classifiers[-1][1][1])
            except IndexError as e:
                max_events_version.append(pd.Series([], name='CorrectProba'))
        max_events[version_prefix_x] = max_events_version
    MetaData().np_saver(max_events, 'files/max_events.pkl')

for version_prefix_x, version_prefix_y in [('GatingPreBool_Stimulus_WindowGatingClassify_PFC', 'Stimulus_StageGatingCenteredMemoryPostDist1Only_WindowInterferenceClassify_PFC'),
                                           ('GatingPreBool_Stimulus_WindowGatingClassify_Stri', 'Stimulus_StageGatingCenteredMemoryPostDist1Only_WindowInterferenceClassify_PFC'),
                                           ('GatingPreBool_Stimulus_WindowGatingClassify_PFC', 'GatedStimulus_StageGatingCenteredSensoryGatingOnly_WindowMemoryClassify_PFC'),
                                           ('GatingPreBool_Stimulus_WindowGatingClassify_PFC', 'Stimulus_StageGatingCenteredMemoryPostDist1Only_WindowInterferenceClassify_PFC')]:
    # fig, ax = plt.subplots(6, 7, figsize=(12, 9))
    # [axi.set_axis_off() for axi in ax.ravel()]
    for sess_ind in range(42):

        print(sess_ind)
        try:

            events_classifier_x = max_events[version_prefix_x][sess_ind]
            events_classifier_x.index.set_names(MetaData().preproc_imports['events']['index'], inplace=True)
            events_classifier_x.rename(version_prefix_x, inplace=True)
            x = db.tables['events_conditions'].merge(events_classifier_x, left_index=True, right_index=True)[['GatingCondSpecialized', version_prefix_x]]
            x = x.loc[x['GatingCondSpecialized'].eq('Gating')]

            events_classifier_y = max_events[version_prefix_y][sess_ind]
            events_classifier_y.index.set_names(MetaData().preproc_imports['events']['index'], inplace=True)
            events_classifier_y.rename(version_prefix_y, inplace=True)
            y = db.tables['events_conditions'].merge(events_classifier_y, left_index=True, right_index=True)[[version_prefix_y]]

            xy = pd.merge(x.reset_index(), y.reset_index(), on=MetaData().preproc_imports['trials']['index'])

            ax_i = ax.ravel()[sess_ind]
            ax_i.scatter(xy[version_prefix_x], xy[version_prefix_y])
            ax_i.axes.xaxis.set_ticklabels([])
            ax_i.axes.yaxis.set_ticklabels([])
            xlim = ax_i.get_xlim()
            ylim = ax_i.get_ylim()
            ax_i.vlines(0.5, ylim[0], ylim[1])
            ax_i.hlines(0.25, xlim[0], xlim[1])

        except ValueError as e:
            pass

    if (version_prefix_x, version_prefix_y)  == ('GatingPreBool_Stimulus_WindowGatingClassify_PFC', 'Stimulus_StageGatingCenteredMemoryPostDist1Only_WindowInterferenceClassify_PFC'):
        fig.suptitle('Gating PFC ~ Gating +1 Late Memory')
    elif (version_prefix_x, version_prefix_y)  == ('GatingPreBool_Stimulus_WindowGatingClassify_Stri', 'Stimulus_StageGatingCenteredMemoryPostDist1Only_WindowInterferenceClassify_PFC'):
        fig.suptitle('Gating Stri ~ Gating +1 Late Memory')
    elif (version_prefix_x, version_prefix_y)  == ('GatingPreBool_Stimulus_WindowGatingClassify_PFC', 'GatedStimulus_StageGatingCenteredSensoryGatingOnly_WindowMemoryClassify_PFC'):
        fig.suptitle('Gating PFC ~ Gating Late Memory')
    elif (version_prefix_x, version_prefix_y)  == ('GatingPreBool_Stimulus_WindowGatingClassify_PFC', 'Stimulus_StageGatingCenteredMemoryPostDist1Only_WindowInterferenceClassify_PFC'):
        fig.suptitle('Gating PFC ~ Gating +1 Interference Sensory')


# Gonzo 6/15, .75, .5/.6
# Oscar 9/15, .75, .6


### classification accuracy

from rec_utils import *
import matplotlib.pyplot as plt
from scipy.stats import binom_test
from scipy.stats import pearsonr
from scipy.stats import spearmanr


db = DataBase(['sessions', 'units', 'events', 'conditions'])
sessions = db.tables['sessions']
md = MetaData()
version_str_list = ['GatingPreBool_Stimulus_WindowGatingClassify_PFC',
                    'GatingPreBool_Stimulus_WindowGatingClassify_Stri',
                    'GatedStimulus_StageGatingCenteredSensoryGatingOnly_WindowInterferenceClassify_PFC',
                    'GatedStimulus_StageGatingCenteredSensoryGatingOnly_WindowMemoryClassify_PFC',
                    'Stimulus_StageGatingCenteredMemoryPostDist1Only_WindowInterferenceClassify_PFC',
                    'GatedStimulus_StageGatingCenteredSensoryPostDist1Only_WindowInterferenceClassify_PFC',
                    'Stimulus_StageGatingCenteredMemoryPostDist1Only_WindowMemoryClassify_PFC',
                    'GatedStimulus_StageGatingCenteredSensoryPostDist1Only_WindowMemoryClassify_PFC']
key = ('15', '0.75', '0.6', 'none', 'multiple', 'lda')

result_dict_matrix = np.empty((8, 42), dtype=bool)
result_df_dict = {}

for v_ii, version_prefix_x in enumerate(version_str_list):
    df_list = []
    for sess_ii in range(42):
        print(v_ii, sess_ii)
        version_str = version_prefix_x + '_{0:d}'.format(sess_ii)
        try:
            target_filename = md.proc_dest_path(path.join('BehavioralUnits', 'GatingCorrelation', 'consolidate', version_str),
                                                'result_dict.pkl')
            result_dict = md.np_loader(target_filename)
            chance = 0.5 if 'GatingPreBool' in version_prefix_x else 0.25
            alpha = 0.05
            p = binom_test(result_dict[key][1].gt(chance).sum(), result_dict[key][1].gt(chance).count(), alternative='greater')

            result_dict_matrix[v_ii, sess_ii] = p < alpha
            if result_dict_matrix[v_ii, sess_ii]:
                df_list.append(result_dict[key][1])
        except KeyError as e:
            result_dict_matrix[v_ii, sess_ii] = False
    result_df_dict[version_prefix_x] = pd.concat(df_list)

md.np_saver((result_dict_matrix, result_df_dict, version_str_list), 'files/gating_correlation_sign_session.pkl')


# plotting
(result_dict_matrix, result_df_dict, version_str_list) = md.np_loader('files/gating_correlation_sign_session_full.pkl')
subject = 'Oscar'
fig, ax = plt.subplots(3, 3, figsize=(12, 9))
[axi.set_axis_off() for axi in ax.ravel()]
for v_ii, (version_prefix_x, version_prefix_y) in enumerate(
        [('GatingPreBool_Stimulus_WindowGatingClassify_PFC', 'GatedStimulus_StageGatingCenteredSensoryGatingOnly_WindowInterferenceClassify_PFC'),
         ('GatingPreBool_Stimulus_WindowGatingClassify_PFC', 'GatedStimulus_StageGatingCenteredSensoryGatingOnly_WindowMemoryClassify_PFC'), #,
         ('GatingPreBool_Stimulus_WindowGatingClassify_PFC', 'Stimulus_StageGatingCenteredMemoryPostDist1Only_WindowInterferenceClassify_PFC'),
         ('GatingPreBool_Stimulus_WindowGatingClassify_Stri', 'GatedStimulus_StageGatingCenteredSensoryGatingOnly_WindowInterferenceClassify_PFC'),
         ('GatingPreBool_Stimulus_WindowGatingClassify_Stri', 'GatedStimulus_StageGatingCenteredSensoryGatingOnly_WindowMemoryClassify_PFC'),
         ('GatingPreBool_Stimulus_WindowGatingClassify_Stri', 'Stimulus_StageGatingCenteredMemoryPostDist1Only_WindowInterferenceClassify_PFC'),
         ('GatingPreBool_Stimulus_WindowGatingClassify_PFC', 'GatingPreBool_Stimulus_WindowGatingClassify_Stri')]):

    try:

        events_classifier_x = result_df_dict[version_prefix_x]
        events_classifier_x.index.set_names(md.preproc_imports['events']['index'], inplace=True)
        events_classifier_x.rename(version_prefix_x, inplace=True)
        x = db.tables['events_conditions'].merge(events_classifier_x, left_index=True, right_index=True)[['GatingCondSpecialized', version_prefix_x]]
        gating_x = 'GatingPreBool' in events_classifier_x
        if gating_x:
            x = x.loc[x['GatingCondSpecialized'].eq('Gating')]

        events_classifier_y = result_df_dict[version_prefix_y]
        events_classifier_y.index.set_names(md.preproc_imports['events']['index'], inplace=True)
        events_classifier_y.rename(version_prefix_y, inplace=True)
        y = db.tables['events_conditions'].merge(events_classifier_y, left_index=True, right_index=True)[['GatingCondSpecialized', version_prefix_y]]
        gating_y = 'GatingPreBool' in events_classifier_y
        if gating_y:
            y = y.loc[y['GatingCondSpecialized'].eq('Gating')]

        xy = pd.merge(x.reset_index(), y.reset_index(), on=MetaData().preproc_imports['trials']['index'])
        xy = xy.loc[xy.reset_index()['Session'].apply(lambda x: subject in x)]
        xy[version_prefix_x] = xy[version_prefix_x].apply(np.exp)
        xy[version_prefix_y] = xy[version_prefix_y].apply(np.exp)
        means = xy[['Session', version_prefix_x, version_prefix_y]].groupby('Session').transform('mean')
        xy[version_prefix_x] = (xy[version_prefix_x] - means[version_prefix_x])
        xy[version_prefix_y] = (xy[version_prefix_y] - means[version_prefix_y])

        # test
        print(pearsonr(xy[version_prefix_x], xy[version_prefix_y]), spearmanr(xy[version_prefix_x], xy[version_prefix_y]))

        # visualization
        xy = xy.sample(min(2000, len(xy)))
        marker = xy['Session'].apply(lambda x: 'o' if 'Oscar' in x else 'v')
        color = xy['Session'].astype('category').cat.codes

        ax_i = ax.ravel()[v_ii]
        ax_i.scatter(xy[version_prefix_x], xy[version_prefix_y], c=color)
        ax_i.axes.xaxis.set_ticklabels([])
        ax_i.axes.yaxis.set_ticklabels([])
        xlim = ax_i.get_xlim()
        ylim = ax_i.get_ylim()
        # ax_i.vlines(0.5 if gating_x else 0.25, ylim[0], ylim[1])
        # ax_i.hlines(0.5 if gating_y else 0.25, xlim[0], xlim[1])
        ax_i.vlines(0, ylim[0], ylim[1], color='#520d18')
        ax_i.hlines(0, xlim[0], xlim[1], color='#520d18')
        if 'PFC' in version_prefix_x: x_str = 'PFC'
        elif 'Stri' in version_prefix_x: x_str = 'Stri'

        if 'GatingOnly' in version_prefix_y: y_str1 = 'Gated'
        elif 'PostDist1Only' in version_prefix_y: y_str1 = 'PostDist'
        elif 'GatingPreBool' in version_prefix_y and 'Stri' in version_prefix_y: y_str1 = 'Stri'
        else: y_str1 = ''

        if 'Interference' in version_prefix_y: y_str2 = 'Early'
        elif 'Memory' in version_prefix_y: y_str2 = 'Late'
        else: y_str2 = ''

        ax_i.set_title('_'.join([x_str, y_str1 + y_str2]))

    except ValueError as e:
        pass

# RESULTS                   FILTERED                                UNFILTERED
#                           GONZO           OSCAR                   GONZO           OSCAR
# PFC ~ Str[]                                                       *** **          *** **
# PFC ~ Gated[300-600]	    *   ***			***  **                 *** ***         *** ***
# PFC ~ Gated[600-900]		    *   - (NEG)		     -                                  **  ***
# Stri ~ Gated[300-600]					    -    - (NEG)                                * (NEG)
# Stri ~ Gated[600-900]
# PFC ~ Memory[300-600]
# PFC ~ Memory[600-900]
#
# GONZO
# FILTERED
# (-0.024759116803701575, 0.19584098984116635) SpearmanrResult(correlation=-0.02934778211739802, pvalue=0.12519875208680398)
# (0.04448719421677488, 0.010674504650203993) SpearmanrResult(correlation=0.07191230442375748, pvalue=3.6206147621576643e-05)
# (-0.029058976273450596, 0.046948352539426384) SpearmanrResult(correlation=-0.02647414544397174, pvalue=0.07029956645467483)
# (-0.02202574663476638, 0.5871670633636724) SpearmanrResult(correlation=-0.007319655598143128, pvalue=0.8568277453988096)
# (-0.021401577426200837, 0.4639934824157921) SpearmanrResult(correlation=-0.008775083210484202, pvalue=0.764006338351698)
# (0.1423157618995755, 0.13622114453672865) SpearmanrResult(correlation=0.14565312647746556, pvalue=0.1271777219172474)
# (-0.06953080658566976, 0.4683672256126994) SpearmanrResult(correlation=-0.060881491327527106, pvalue=0.5255885090441813)
# UNFILTERED
# (-0.030788646688738326, 0.0007265603776737713) SpearmanrResult(correlation=-0.02653434237422616, pvalue=0.0035886276848411977)
# (0.0808910294381627, 2.847993384695422e-16) SpearmanrResult(correlation=0.05907797443033464, pvalue=2.3825868582661904e-09)
# (-0.000515964839781141, 0.958458430381815) SpearmanrResult(correlation=0.010422567627104149, pvalue=0.29269971575669473)
# (-0.009060877443250251, 0.5065140280812959) SpearmanrResult(correlation=-0.0027889236225403747, pvalue=0.8379947300355755)
# (-0.019415450016806057, 0.1545913008429661) SpearmanrResult(correlation=-0.01753233624013744, pvalue=0.1986485008325501)
# (0.03533903938381378, 0.4034842715212048) SpearmanrResult(correlation=0.044461121374538524, pvalue=0.29314276710015813)
# (-0.06608987866359008, 0.22349849065840616) SpearmanrResult(correlation=-0.034974031388499395, pvalue=0.5197943437884652)
#
#
#
# OSCAR
# FILTERED
# (0.019994769187952226, 0.21106672848163943) SpearmanrResult(correlation=0.024839294309414936, pvalue=0.1202468033103246)
# (0.09758191262940899, 2.995903986997309e-11) SpearmanrResult(correlation=0.041119986151548306, pvalue=0.005183949189383958)    ***     **
# (0.014442759135441598, 0.27662889072170693) SpearmanrResult(correlation=0.02557278688639657, pvalue=0.054038872832348664)                -
# (-0.06959249162748049, 0.07488010190777213) SpearmanrResult(correlation=-0.06850795392578662, pvalue=0.07953787832523887)         -       -   NEG
# (0.0005726537571820553, 0.9809343057412016) SpearmanrResult(correlation=0.03655315511999385, pvalue=0.12702978120034378)
# (-0.04149242143682101, 0.6366548333032161) SpearmanrResult(correlation=-0.04949732869896944, pvalue=0.5730142607223523)
# (-0.06734071003096842, 0.4429634460390993) SpearmanrResult(correlation=-0.072842151456138, pvalue=0.40651834656455177)
# UNFILTERED
# (0.04739393369012355, 1.9275981629366866e-06) SpearmanrResult(correlation=0.026723726075711947, pvalue=0.0072866114601372305)
# (0.10900853719719426, 2.2093416519974994e-19) SpearmanrResult(correlation=0.055676222596044614, pvalue=4.4791609731781954e-06)
# (0.03199135374989238, 0.008419444852055182) SpearmanrResult(correlation=0.046077559424773144, pvalue=0.0001470866989991429)
# (0.0038119349353617333, 0.7901094503978382) SpearmanrResult(correlation=-0.02946883587566269, pvalue=0.03958034796431191)
# (0.0012380469552243022, 0.9311115706412225) SpearmanrResult(correlation=-0.012734478388486908, pvalue=0.37388581481618677)
# (0.03733254224661472, 0.5398219216217202) SpearmanrResult(correlation=0.03637067327894531, pvalue=0.5503233247800625)
# (-0.041355706904762125, 0.5836223724174814) SpearmanrResult(correlation=-0.05364888106974614, pvalue=0.4769386121425724)

from rec_utils import *
import matplotlib.pyplot as plt
from scipy.stats import binom_test
from scipy.stats import pearsonr
from scipy.stats import spearmanr

# plotting
# fig, ax = plt.subplots(3, 3, figsize=(12, 9))
# [axi.set_axis_off() for axi in ax.ravel()]

md = MetaData()
db = DataBase(['sessions', 'units', 'events', 'conditions'])
events_conditions = db.tables['events_conditions']
# create arg version lists
args_counts_thr = ['counts_thr=12']
args_sess_ratio = ['sess_ratio={0:s}'.format(ratio) for ratio in ['0.75']]
args_units_ratio = ['units_ratio={0:s}'.format(ratio) for ratio in ['0.6']]
args_scaler = ['scaler=none']
args_imputer = ['imputer=multiple']
args_classifier = ['classifier=lda']

args_version_list_x = []
for class_arg, balance, area_list in [('GatingPreBool', 'Stimulus', 'PFC'),
                                      ('GatingPreBool', 'Stimulus', 'Stri')]:
    args_class = ['class={0:s}'.format(class_arg)]
    args_balance = ['balance={0:s}'.format(balance)]
    args_area_list = ['area_list={0:s}'.format(area_list)]

    args_version_list_x.extend(list(map(list, list(product(args_class, args_balance, args_counts_thr, args_area_list,
                                                         args_sess_ratio, args_units_ratio, args_scaler,
                                                         args_imputer, args_classifier)))))
args_version_list_y = []
for class_arg, balance, area_list in [('GatedStimulus', 'StageGatingCenteredSensoryGatingOnly', 'PFC'),
                                      ('GatedStimulus', 'StageGatingCenteredSensoryPostDist1Only', 'PFC'),
                                      ('Stimulus', 'StageGatingCenteredMemoryPostDist1Only', 'PFC')]:
    args_class = ['class={0:s}'.format(class_arg)]
    args_balance = ['balance={0:s}'.format(balance)]
    args_area_list = ['area_list={0:s}'.format(area_list)]
    args_version_list_y.extend(list(map(list, list(product(args_class, args_balance, args_counts_thr, args_area_list,
                                                           args_sess_ratio, args_units_ratio, args_scaler,
                                                           args_imputer, args_classifier)))))

full_results = {}
# for all combinations of arg

for args_version_x, args_version_y in product(args_version_list_x, args_version_list_y):

    print(args_version_x, args_version_y)
    def process_results(args_fr, args_version):
        events_series_list = []
        valid_list = []
        for args_subject in [str(sess) for sess in range(42)]:
            try:
                version = job_scheduler(args_version)
                version['fr'] = args_fr
                version['subject'] = args_subject
                classifier = ClassificationAnalysis(DataBase(['units']), version)
                target_filename = classifier.get_path_base('events', classifier.get_train_test_session_stem())
                events_series = md.np_loader(target_filename)

                chance = 0.5
                alpha = 0.05
                p = binom_test(np.exp(events_series).gt(chance).sum(), np.exp(events_series).gt(chance).count(),
                               alternative='greater')
                valid = p < alpha
                valid_list.append(valid)
                events_series_list.append(events_series)
            except:
                pass
        return (valid_list, events_series_list)

    x = {}
    for args_fr in ['200_400']:
        print(args_fr)
        x[args_fr] = process_results(args_fr, args_version_x)

    y = {}
    for args_fr in ['200_400', '400_600', '600_800', '800_1000']:
        print(args_fr)
        y[args_fr] = process_results(args_fr, args_version_y)

    for args_fr_x, args_fr_y in product(['200_400'], ['200_400', '400_600', '600_800', '800_1000']):

        try:
            for es in x[args_fr_x][1]:
                es.index.set_names(md.preproc_imports['events']['index'], inplace=True)
                es.rename('x', inplace=True)
                es = events_conditions.merge(es, left_index=True, right_index=True)[['GatingCondSpecialized', 'x']]
                gating_x = 'GatingPreBool' in job_scheduler(args_version_x)['class']
                if gating_x:
                    es = es.loc[es['GatingCondSpecialized'].eq('Gating')]
            x_df = pd.concat([es - es.mean() for v, es in list(zip(*x[args_fr_x]))], axis=0)

            for es in y[args_fr_y][1]:
                es.index.set_names(md.preproc_imports['events']['index'], inplace=True)
                es.rename('y', inplace=True)
                es = events_conditions.merge(es, left_index=True, right_index=True)[['GatingCondSpecialized', 'y']]
                gating_y = 'GatingPreBool' in job_scheduler(args_version_y)['class']
                if gating_y:
                    es = es.loc[es['GatingCondSpecialized'].eq('Gating')]
            y_df = pd.concat([es - es.mean() for v, es in list(zip(*y[args_fr_y]))], axis=0)

            xy = pd.merge(x_df.reset_index(), y_df.reset_index(), on=MetaData().preproc_imports['trials']['index'])
            print(pearsonr(xy['x'], xy['y']),
                  spearmanr(xy['x'], xy['y']))
            # take pearson/spearman correlation

            full_results[('_'.join(args_version_x), '_'.join(args_version_y))] = (xy, pearsonr(xy['x'], xy['y']), spearmanr(xy['x'], xy['y']))
        except:
            full_results[('_'.join(args_version_x), '_'.join(args_version_y))] = np.nan


# 'PFC', 'Gated'
# (0.02555999499330899, 0.00046038712586561885)
# (0.008212342881730611, 0.2604726122687524)
# (-0.012587483487292942, 0.08456240168948209)
# (0.0021814870420572784, 0.7650103449903518)

# 'PFC', 'Memory'
# (0.0076357392690432, 0.708432707408085)
# (-0.032727321662745605, 0.10888319437985491)
# (0.002869657216844411, 0.8882331374305634)
# (-0.003960794329461592, 0.8461926242088633)

# 'PFC', 'Sensory'
# (-0.030947089076851352, 0.13162064622900146)
# (-0.007303809064163956, 0.7220195693464615)
# (-0.019333303194306505, 0.3463049493523551)
# (0.0036452395086926867, 0.8590731997389619)

# 'Stri', 'Gated'
# (-0.012756832551151742, 0.17901374256066596)
# (0.022830355748902197, 0.016165692842998593)
# (-0.0006823432195703453, 0.9427015816337526)
# (-0.010270800030206985, 0.2792946998873958)

# 'Stri', 'Memory'
# (-0.07412106714781987, 0.0058554158321718975) **
# (-0.005550190643472678, 0.8367362789956214)
# (0.021387024444274174, 0.42710715490400863)
# (0.02955511172485504, 0.27239378166997313)

# 'Stri', 'Sensory'
# (0.0028987223020363567, 0.9148840433593582)
# (-0.00666075082426955, 0.8059965900815073)
# (7.460616779907167e-05, 0.9978051549416693)
# (0.0006155711938930647, 0.981891990024923)

from statsmodels.stats.multitest import multipletests
multipletests([.9818, .9978, .8059, .9148, .2723, .4271, .8367, .0058, .2792, .9427, .0161, .1790, \
               .8590, .3463, .7220, .1316, .8461, .8882, .1088, .7084, .7650, .0845, .2604, .0004])