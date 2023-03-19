# PostDistractor

import seaborn as sns
from scipy.stats import ttest_rel
from rec_utils import *
version_aov = 'GatedStimulusPostDistractorCentered'
version_fr = 'WindowDelayLate'
v_aov_params = anova_version_aov_params(version_aov, version_fr)
x_factors = v_aov_params['x_factors']io
selection_dict = v_aov_params['selection_dict']
v_fr_params = version_fr_params(version_fr)
t_start = v_fr_params['t_start']
t_end = v_fr_params['t_end']
timebin = v_fr_params['timebin']
timestep = v_fr_params['timestep']

# load data and init vars, tables, and slices
md = MetaData()
db = md.db_base_loader(['physiology'])
physiology = db['physiology']
target_filename = md.proc_dest_path(path.join('BehavioralUnits', 'Anova', version_aov,
                                              behunit_params_str(version_fr, timebin, timestep, t_start, t_end),
                                              'summarize'),
                                    'physiology_dict.pkl')
physiology_dict = md.np_loader(target_filename)
db = DataBase()
units = db.tables['units']



temp = [np.array([sel_ph_val[x_factors[0]]['zscores'][0]
                  for u_ind, sel_ph_val
                  in sel_physiology.items()
                  if units.loc[u_ind].Area == 'PFC'
                  and units.loc[u_ind].UnitNum != 0
                  and units.loc[u_ind].RatingCode != 7
                  and bool(sel_ph_val[x_factors[0]]['zscores'])])
        for selection, sel_physiology
        in physiology_dict.items()]
print(ttest_rel(temp[0], temp[1]))
print(ttest_rel(temp[0], temp[2]))
print(ttest_rel(temp[0], temp[3]))
print(ttest_rel(temp[1], temp[2]))
print(ttest_rel(temp[1], temp[3]))
print(ttest_rel(temp[2], temp[3]))
palette = ['#dbb356', '#3c3a7a', '#c9834d', '#a1673a']
data = pd.DataFrame(np.array(temp).transpose(), columns=['Gating-1', 'Gating', 'Gating+1', 'Gating+2'])
# ax = sns.boxplot(data=temp, palette=['#dbb356', '#3c3a7a', '#c9834d', '#a1673a'], saturation=1)
ax = sns.catplot(data=data, palette=palette, saturation=1, kind='boxen')
# ax = sns.catplot(data=data.agg(lambda x: x - data['Gating-1']), palette=palette, saturation=1, kind='boxen')
ax = ax.ax_i
# ax.set_xticklabels(['Gating-1', 'Gating', 'Gating+1', 'Gating+2'])
ax.set_ylabel('z-scored $\omega^2$')
ax.set_position([0.2, 0.1, 0.75, 0.8])

####################
####################
####################
from rec_utils import *
from aov_stats import *
version_aov = 'GeneralizedGating'
import seaborn as sns
import matplotlib.pyplot as plt

# DoSSI scatters and histograms
for version_fr in ['WindowSampleShift', 'WindowDelayShift']:
    if version_fr == 'WindowDelayShift':
        stim_str_levels = ['S11D_S12D_S21D_S22D', 'S11D', 'S12D', 'S21D', 'S22D']

    md = MetaData()
    db = md.db_base_loader(['units', 'physiology'])
    units, physiology = db['units'], db['physiology']

    v_fr_params = anova_version_fr_params(version_fr)
    t_start = v_fr_params['t_start']
    t_end = v_fr_params['t_end']
    timebin = v_fr_params['timebin']
    timestep = v_fr_params['timestep']
    target_filename = md.proc_dest_path(path.join('BehavioralUnits', 'Anova', version_aov,
                                                  behunit_params_str(version_fr, timebin, timestep, t_start, t_end),
                                                  'selectivity_index'),
                                        'physiology_df.pkl')
    physiology_df = md.np_loader(target_filename)

    area_df = []
    f, ax = plt.subplots(1, 3)
    f.set_size_inches(12, 4)
    for area_i, area in enumerate(['PFC', 'Stri', 'IT']):
        area_units = units.loc[units['Area'].eq(area) & units['UnitNum'].ne(0) & units['RatingCode'].ne(7)]
        area_physiology = physiology_df.loc[area_units.index]
        df = area_physiology[['DepthOfSelectivityIndex_Distractor', 'DepthOfSelectivityIndex_Gating']]
        df = df.rename({'DepthOfSelectivityIndex_Distractor': 'Dist', 'DepthOfSelectivityIndex_Gating': 'Gating'}, axis=1)
        df['Area'] = area
        area_df.append(df)

        sns.regplot(data=df, x='Dist', y='Gating', color='#a61629', ax=ax[area_i])
        ax[area_i].set_title(area)
        ax[area_i].plot([0, 1], [0, 1], c='#2B2B2B')

    df = pd.concat(area_df, axis=0)
    df_melt = df.reset_index().melt(id_vars=['Session', 'ChanNum', 'UnitNum', 'Area'], var_name='Period', value_name='DoSSI')
    plt.figure()
    sns.violinplot(data=df_melt, x='Area', hue='Period', y='DoSSI', split=True, palette=['#dbb356', '#3c3a7a'])



########
########

# anova testing for nested anova test

ratings = np.array([9,8,6,8,10,4,6,5,7,7, 7,9,6,6,6,11,6,3,8,7, 11,13,8,6,14,11,13,13,10,11, 12,11,16,11,9,23,12,10,19,11, 10,19,14,5,10,11,14,15,11,11,
                    8,6,4,6,7,6,5,7,9,7, 10,7,8,10,4,7,10,6,7,7, 14,11,18,14,13,22,17,16,12,11, 20,16,16,15,18,16,20,22,14,19, 21,19,17,15,22,16,22,22,18,21])
ratings = ratings.reshape(1, -1)
therapist = np.tile(np.arange(10), (int(100 / 10), 1)).transpose().reshape(1, -1)
gender = np.tile(np.arange(2), (int(100 / 2), 1)).transpose().reshape(1, -1)
customer = np.tile(np.arange(10), (int(100 / 10), 1)).reshape(1, -1)
df = pd.DataFrame(np.concatenate((ratings, therapist, gender, customer)).transpose(), columns=['Rating', 'Therapist', 'Gender', 'Customer'])


##############################
##############################

# Unit counts, Venn diagrams

from rec_utils import *
from matplotlib_venn import venn3

db = DataBase()
units = db.tables['units']
anova_sens = Anova(DataBase([]), {'aov': 'PresentedStimulus', 'fr': 'ConcatFactor'})
anova_mem = Anova(DataBase([]), {'aov': 'GatedStimulus', 'fr': 'ConcatFactor'})
x_mem = 'RuleStimCategory'
x_sens = 'StageStimSpecialized'

units_slice = units.loc[units['Area'].eq('PFC') & units['UnitNum'].ne(0) & units['RatingCode'].ne(7)]
mem_units_inds = anova_mem.get_all_units('Gating', x_mem)
sens_units_inds = anova_sens.get_all_units('Gating', x_sens)

gating_mem_units_inds = anova_mem.get_selective_units('Gating', x_mem)
gating_sens_units_inds = anova_sens.get_selective_units('Gating', x_sens)
post_mem_units_inds = anova_mem.get_selective_units('PostDist', x_mem)
post_sens_units_inds = anova_sens.get_selective_units('PostDist', x_sens)
pre_mem_units_inds = anova_mem.get_selective_units('PreDist', x_mem)
pre_sens_units_inds = anova_sens.get_selective_units('PreDist', x_sens)

mem_units = set(units_slice.index).intersection(mem_units_inds)
sens_units = set(units_slice.index).intersection(sens_units_inds)
gating_mem_units = set(units_slice.index).intersection(gating_mem_units_inds)
gating_sens_units = set(units_slice.index).intersection(gating_sens_units_inds)
post_mem_units = set(units_slice.index).intersection(post_mem_units_inds)
post_sens_units = set(units_slice.index).intersection(post_sens_units_inds)
pre_mem_units = set(units_slice.index).intersection(pre_mem_units_inds)
pre_sens_units = set(units_slice.index).intersection(pre_sens_units_inds)

a = post_mem_units
b = post_sens_units
c = gating_sens_units

abc = set(a).intersection(b).intersection(c)
ab = set(a).intersection(b).difference(abc)
ac = set(a).intersection(c).difference(abc)
bc = set(b).intersection(c).difference(abc)
ax = set(a).difference(abc).difference(ab).difference(ac)
bx = set(b).difference(abc).difference(ab).difference(bc)
cx = set(c).difference(abc).difference(ac).difference(bc)

for _ in range(4):
    venn3(list(map(len, [ax, bx, ab, cx, ac, bc, abc])),
          set_labels=['PostDist Memory', 'PostDist Sensory', 'Gated Stimulus'],
          set_colors=['#e08e4f', '#e06132', '#3c3a7a'])


##############################


#Presented stimulus count significance
# PFC, IT, Stri x PreDist, Gating, PostDist, All, Any

from scipy.stats import binom_test
from statsmodels.stats.multitest import multipletests
multipletests([binom_test(x, n, 0.05) for x, n, p in zip([104, 162, 131, 31, 276, 18, 14, 18, 8, 29, 20, 37, 26, 3, 67],
                                                         [897, 897, 897, 897, 897, 102, 102, 102, 102, 102, 394, 394, 394, 394, 394],
                                                         [.05, .05, .05, .000125, .142625, .05, .05, .05, .000125, .142625, .05, .05, .05, .000125, .142625])])
