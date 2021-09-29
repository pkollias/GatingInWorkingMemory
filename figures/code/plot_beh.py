from rec import *

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# from sklearn.linear_model import LogisticRegression
# import statsmodels.api as sm


# init params
plots = {}
fsize = (6, 4)
barWidth = 0.8
hold_color = '#BDC3C7'
release_color = '#1B2631'
correct_color = '#2ECC71'
incorrect_color = '#E74C3C'
ordinal = lambda n: "%d%s" % (n, "tsnrhtdd"[(n / 10 % 10 != 1) * (n % 10 < 4) * n % 10::4])
def repeat_responses(row):
    if not row['GatingCondExtended'] in ['PreDist', 'PostDist', 'Target']:
        return np.nan
    elif row['GatingCondExtended'] == 'Target':
        return 'Target'
    elif row['GatingCondExtended'] == 'Gating':
        return 'Gating'
    else:
        return '{0:s}'.format(ordinal(int(row['StimOccurrence'])))



# data collection
md = MetaData()
db = md.db_base_loader(['sessions', 'trials', 'events', 'conditions'])
sessions, trials, events, conditions = db['sessions'], db['trials'], db['events'], db['conditions']
trials_index = md.preproc_imports['trials']['index']
events_index = md.preproc_imports['events']['index']
events_conditions = pd.merge(events.reset_index(drop=True), conditions, on=events_index)

# data wrangling
trials['Correct'] = trials['StopCondition'].eq(1)
trials['BreakFix'] = trials['StopCondition'].eq(-3)
stage_index_pair = events_conditions.reset_index()['StageIndex'].apply(lambda x: np.ceil(x / 2))
events_conditions['StageIndexPair'] = stage_index_pair
events_conditions.set_index(events_index, inplace=True)

# trials_last_event = events_conditions.groupby(trials_index).tail(1).set_index(trials_index, drop=False)


for subject in ['Oscar', 'Gonzo']:

    trials_subject = trials.loc[sessions.loc[sessions['Subject'].eq(subject)].index]
    events_conditions_subject = events_conditions.loc[sessions.loc[sessions['Subject'].eq(subject)].index]

    valid = trials_subject.loc[trials['StopCondition'].isin([1, -3, -5, -6, -7, -8])]
    attempted = trials_subject.loc[trials['StopCondition'].isin([1, -5, -6, -7, -8])]
    incorrect = attempted.loc[~attempted['Correct']]
    early = incorrect.loc[incorrect['StopCondition'].eq(-7)]

    # compress events_conditions
    events_compressed = events_conditions_subject.reset_index().drop_duplicates(
        subset=trials_index + ['StageIndexPair'], keep='last').set_index(trials_index)
    attempted_compressed = events_compressed.loc[attempted.index]
    attempted_compressed['GatingOccurrenceCategory'] = attempted_compressed.apply(repeat_responses, axis=1)

    # intermediate results
    sc_series = (attempted.StopCondition.value_counts() / len(attempted)).loc[[1, -7, -6]]

    ax, fig = {}, {}


    # Response Types
    plot_label = 'AttemptedOutcome'

    c
    data = sc_series.to_list()
    colors = ['#14b353', '#ad2415', '#919191']
    labels = ['CORRECT', 'EARLY', 'LATE']
    explode = [0.05 for _ in range(3)]
    wedges, texts = plt.pie(data, colors=colors, labels=['', '', ''], explode=explode, wedgeprops=dict(width=1),
                            startangle=90 + data[-1] * 360)

    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    kw = dict(arrowprops=dict(arrowstyle="-"),
              bbox=bbox_props, zorder=0, va="center")

    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1) / 2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        ax.annotate('{0:.2f}%'.format(100 * data[i]), xy=(x, y), xytext=(1.35 * np.sign(x), 1.4 * y),
                    horizontalalignment=horizontalalignment, **kw)
    if subject == 'Gonzo':
        plt.legend(wedges, labels, frameon=False)
        ax.get_legend().set_bbox_to_anchor((1, .75))

    ax.get_figure().savefig('figures/Behavior/StopCondition_{0:s}.png'.format(subject), format='png')
    # ax.set_position([0.3, 0.3, 0.4, 0.4])
    # fig[plot_label] = plt.figure(figsize=(7, 7))
    # ax[plot_label] = sc_series.plot.bar()
    # ax[plot_label].set_xticklabels(('CORRECT', 'EARLY', 'LATE'))
    # ax[plot_label].set_ylim((0, 1))
    # ax[plot_label].set_ylabel('Ratio over Attempted Trials')
    # ax[plot_label].set_title('Attempted Trials Outcome Distribution')


    # accuracy by number of distractors
    plot_label = 'DistAccuracy'
    fig[plot_label] = plt.figure(figsize=fsize)
    nn = attempted.groupby(['NumPre', 'NumPost']).size().rename('nn')
    pp = attempted.groupby(['NumPre', 'NumPost']).mean()['Correct'].rename('pp')
    npp = (1 - attempted.groupby(['NumPre', 'NumPost']).mean()['Correct']).rename('npp')
    dist_acc = pd.concat([nn, pp, npp], axis = 1)
    dist_acc['std'] = dist_acc.apply(lambda row: np.sqrt(row['pp'] * row['npp'] / row['nn']), axis=1)
    labels = ['{0:.2f}{1:s}{2:.2f}\n(N={3:d})'.format(vals['pp'], u"\u00B1", 2 * vals['std'], int(vals['nn'])) for _, vals
              in dist_acc.iterrows()]
    data = attempted.pivot_table(values='Correct', index='NumPre', columns='NumPost', aggfunc=np.mean)
    ax[plot_label] = sns.heatmap(
        attempted.pivot_table(values='Correct', index='NumPre', columns='NumPost', aggfunc=np.mean), vmin=0.75, vmax=1,
        annot=np.reshape(labels, data.get_current_shape), fmt='', annot_kws={"fontsize":9})
    ax[plot_label].set_title('Accuracy by Number of Distractors')


    # break fixation by number of distractors
    plot_label = 'DistBreakFix'
    fig[plot_label] = plt.figure(figsize=fsize)
    ax[plot_label] = sns.heatmap(valid.pivot_table(values='BreakFix', index='NumPre', columns='NumPost', aggfunc=np.mean), vmin=0, vmax=1, annot=True, fmt=".2f")
    ax[plot_label].set_title('Break Fixation Rate by Number of Distractors')


    # accuracy by stim/cue
    plot_label = 'StimAccuracy'
    fig[plot_label] = plt.figure(figsize=fsize)
    ax[plot_label] = sns.barplot(x="RuleStimCategory", y="Correct", hue='RuleCueCategory', data=attempted)
    ax[plot_label].set_title('Accuracy by Gated Stimulus Identity')
    ax[plot_label].set_ylim((0.75, 1))
    ax[plot_label].legend(loc='upper right', ncol=1, fontsize='xx-small')


    # bar status by stimulus repetition
    plot_label = 'StimulusRepetition'
    df = attempted_compressed.pivot_table(values='Catch', index='GatingOccurrenceCategory', columns='BarStatus',
                                          aggfunc='count').fillna(0)
    totals = [i + j for i, j in zip(df['HOLD'], df['RELEASE'])]
    release_bars = [i / j * 100 for i, j in zip(df['RELEASE'], totals)]
    hold_bars = [i / j * 100 for i, j in zip(df['HOLD'], totals)]
    names = ['{0:s}\n(N={1:d})'.format(name, int(total)) for name, total in zip(list(df.index), totals)]
    r = list(np.arange(len(names)))
    fig[plot_label], ax[plot_label] = plt.subplots()
    plt.bar(r, release_bars, color=release_color, edgecolor=[incorrect_color for _ in np.arange(len(names) - 1)] + [correct_color], linewidth=1, width=barWidth)
    plt.bar(r, hold_bars, bottom=release_bars, color=hold_color, edgecolor=[correct_color for _ in np.arange(len(names) - 1)] + [incorrect_color], linewidth=1, width=barWidth)
    plt.title('Lever Behavior by Ordinal Occurrence of Stimulus in Sequence')
    plt.xticks(r, names)
    plt.xlabel("Ordinal Occurrence of Stimulus Category")
    plt.ylabel("Stacked Probability of Lever Behavior")
    legend_elements = [Line2D([0], [0], color='#424949', lw=4, label='Release'),
                       Line2D([0], [0], color='#F2F3F4', lw=4, label='Hold'),
                       Patch(facecolor='white', edgecolor='#2ECC71', label='Correct'),
                       Patch(facecolor='white', edgecolor='#E74C3C', label='Incorrect')]
    ax[plot_label].legend(handles=legend_elements, loc='upper left')


    # bar release types
    plot_label = 'BarRelease'
    release_observed = events_compressed.loc[events_compressed['BarStatus'].eq('RELEASE')]
    release_opportunity = attempted_compressed
    gating_stim_values = pd.concat([release_observed['GatingCondStageStimCategory'].value_counts().sort_index().rename('Observed'),
                                   release_opportunity['GatingCondStageStimCategory'].value_counts().sort_index().rename('Opportunity')],
                                  axis=1)
    gating_stim_values.sort_values(by='Opportunity', ascending=False, inplace=True)
    results = gating_stim_values[cols].div(gating_stim_values[cols].sum(axis=0), axis=1)
    labels = list(results.index)
    opportunity_counts = list(results['Opportunity'])
    observed_counts = list(results['Observed'])
    x_vals = np.arange(len(labels)) + 1

    fig[plot_label] = plt.figure(constrained_layout=True)
    spec = fig[plot_label].add_gridspec(ncols=1, nrows=2, width_ratios=[1], height_ratios=[1, 4])
    ax[plot_label] = {1: fig[plot_label].add_subplot(spec[0, 0]),
                      2: fig[plot_label].add_subplot(spec[1, 0])}
    ax1 = ax[plot_label][1]
    ax2 = ax[plot_label][2]
    width = 0.35
    rects11 = ax1.bar(x_vals - width / 2, opportunity_counts, width, label='Opportunity')
    rects12 = ax1.bar(x_vals + width / 2, observed_counts, width, label='Observed')
    rects21 = ax2.bar(x_vals - width / 2, opportunity_counts, width, label='Opportunity')
    rects22 = ax2.bar(x_vals + width / 2, observed_counts, width, label='Observed')
    ax1.set_ylim(.85, 1.)  # outliers only
    ax2.set_ylim(0, .25)  # most of the data
    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax1.xaxis.tick_top()
    ax1.tick_params(labeltop=False)  # don't put tick labels at the top
    ax2.xaxis.tick_bottom()
    d = .015  # how big to make the diagonal lines in axes coordinates
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
    ax2.set_xticklabels([''] + labels)
    plt.xticks(rotation=30)
    ax1.legend()
    ax1.set_title('Observed and Opportunity for Bar Release in Different Event Categories')
    ax2.set_xlabel('Event Category')
    ax2.set_ylabel('Ratio of Events')
    plt.savefig('{0:s}_{1:s}.png'.format(plot_label, subject), format='png')

    for key, axi in ax.items():
        axi.get_figure().savefig('{0:s}_{1:s}.png'.format(key, subject), format='png')



### TASK

# NumRewards by monkey plot
df = pd.merge(trials.reset_index(drop=True), sessions.reset_index(drop=True), on='Session')
sns.boxplot(data=df.loc[df['StopCondition'].eq(1)], x='NumDist', y='NumRewards', hue='Subject')
