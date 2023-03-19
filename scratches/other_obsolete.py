from rec_analyses import *
from dPCA.utils import denoise_mask
from scipy.stats import ttest_ind

areas = ['PFC', 'Stri']
n_areas = 2
n_factors = 10
n_bootstraps = 1000

# ### Bootstrap ### #

args_version_base = ['factor=StimulusGatingPreBool', 'fr=ConcatFactor2', 'counts_thr=20',
'area_list=PFC_Stri', 'subject=Gonzo_Oscar', 'mode=Bootstrap']

fn = 'files/onset_array_bootstrap.pkl'
if path.exists(fn):
    onset_array = pd.read_pickle(fn)
else:
    onset_array = np.empty((n_areas, n_factors, n_bootstraps))
    for area_i, area in enumerate(areas):
        for mode_seed_i in range(n_bootstraps):
            for factor_i in range(n_factors):

                print(area_i, mode_seed_i, factor_i)
                args_version = args_version_base + ['area={0:s}'.format(area)] + ['mode_seed={0:d}'.format(mode_seed_i)]
                version = parse_vars(args_version)
                dpca = DemixedPrincipalComponent(DataBase([]), version)
                target_filename = dpca.get_exec_filename('significance')
                significance = dpca.db.md.np_loader(target_filename)

                significance_timeseries = denoise_mask(significance['g'][factor_i], 4)
                significant_bins = np.where(significance_timeseries)[0]
                num_significant_bins = significant_bins.size
                onset = np.nan if not bool(num_significant_bins) else significant_bins[0]

                onset_array[area_i, factor_i, mode_seed_i] = onset

    with open(fn, 'wb') as f:
        pickle.Pickler(f).dump(onset_array)

df = pd.DataFrame(onset_array.reshape((-1,1)),
                  index=pd.MultiIndex.from_product([range(n_areas), range(n_factors), range(n_bootstraps)],
                                                   names=('Area', 'Factor', 'Bootstrap')))

ax = sns.boxplot(x="Factor", y=0, hue="Area", data=df.reset_index(), palette="muted")
print('PFC:PC1 - Stri:PC1', ttest_ind(df.loc[0, 0, :], df.loc[1, 0, :]))
print('PFC:PC2 - Stri:PC1', ttest_ind(df.loc[0, 1, :], df.loc[1, 0, :]))

# ### Shuffles ### #

args_version_base = ['factor=StimulusGatingPreBool', 'fr=ConcatFactor2', 'counts_thr=20',
'area_list=PFC_Stri', 'subject=Gonzo_Oscar', 'mode=AreaShuffle'] ###

fn = 'files/onset_array_areashuffle.pkl' ###
count_missing = 0
if path.exists(fn):
    onset_array = pd.read_pickle(fn)
else:
    onset_array = np.empty((n_areas, n_factors, n_bootstraps))
    for area_i, area in enumerate(areas):
        for mode_seed_i in range(n_bootstraps):
            for factor_i in range(n_factors):

                print(area_i, mode_seed_i, factor_i)
                args_version = args_version_base + ['area={0:s}'.format(area)] + ['mode_seed={0:d}'.format(mode_seed_i)]
                version = parse_vars(args_version)
                dpca = DemixedPrincipalComponent(DataBase([]), version)
                target_filename = dpca.get_exec_filename('significance')
                try:
                    significance = dpca.db.md.np_loader(target_filename)
                    significance_timeseries = denoise_mask(significance['g'][factor_i], 4)
                    significant_bins = np.where(significance_timeseries)[0]
                    num_significant_bins = significant_bins.size
                    onset = np.nan if not bool(num_significant_bins) else significant_bins[0]
                except:
                    count_missing = count_missing + 1
                    onset = np.nan

                onset_array[area_i, factor_i, mode_seed_i] = onset

    with open(fn, 'wb') as f:
        pickle.Pickler(f).dump(onset_array)

df = pd.DataFrame(onset_array.reshape((-1,1)),
                  index=pd.MultiIndex.from_product([range(n_areas), range(n_factors), range(n_bootstraps)],
                                                   names=('Area', 'Factor', 'Shuffle')))



################################################################################
################################################################################

def evaluate_anova_shuffles(shuffles_series, percentile, shuffles):

    rounding_decimal = 10
    anova_round = lambda x: round(x, rounding_decimal)
    shuffles_index = [shuffle_to_name(ii) for ii in range(shuffles)]

    if len(shuffles_series) == 1:
        valid = False
        if 'observed' in shuffles_series.index:
            omega_sq_observed = shuffles_series.loc['observed']
        else:
            omega_sq_observed = 0
        omega_sq_distr = pd.Series(data=omega_sq_observed, index=shuffles_index)
        threshold = np.inf
        omega_sq_diff = pd.Series(data=-np.inf, index=shuffles_index)
        zscore = 0
    else:
        valid = True
        omega_sq_distr = shuffles_series.apply(anova_round)
        omega_sq_observed = omega_sq_distr.loc[shuffle_to_name(0)]
        shuffle_values = omega_sq_distr.loc[~omega_sq_distr.index.isin([shuffle_to_name(0)])]
        omega_sq_mean = np.mean(shuffle_values)
        omega_sq_std = np.std(shuffle_values)
        threshold = anova_round(np.percentile(shuffle_values, percentile, interpolation='linear'))
        omega_sq_diff = omega_sq_distr.apply(lambda x: anova_round(x - threshold))
        zscore = anova_round((omega_sq_observed - omega_sq_mean) / omega_sq_std) if omega_sq_std > 0 else 0

    return {'valid': valid,
            'omega_sq_distr': omega_sq_distr,
            'omega_sq_observed': omega_sq_observed,
            'threshold': threshold,
            'omega_sq_diff': omega_sq_diff,
            'zscore': zscore}



####################
####################
####################


from rec_analyses import *
from rec import SignalSmoothing

# scatters of gating factors
for area_i, area in enumerate(['PFC', 'Stri']):

    args_version = ['factor=StimulusGatingPreBool', 'fr=ConcatFactor2', 'counts_thr=15',
    'area_list=PFC_Stri', 'subject=Gonzo_Oscar', 'area={0:s}'.format(area), 'mode=Bootstrap', 'mode_seed=0']

    version = parse_vars(args_version)
    dpca = DemixedPrincipalComponent(DataBase([]), version)
    db, md = dpca.db, dpca.db.md
    target_filename = dpca.get_path_base('dpca_obj', dpca.get_exec_stem())
    dpca_obj = md.np_loader(target_filename)

    f, axs = plt.subplots(1, 3)
    f.set_size_inches(16.75, 5)
    for comb_i, comb in enumerate([(0, 1), (0, 2), (1, 2)]):
        ax = axs[comb_i]
        ax.scatter(dpca_obj.D['g'][:, comb[0]], dpca_obj.D['g'][:, comb[1]])
        ax.axis('equal')
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        ax.set_xlabel('PC{0:d}'.format(comb[0] + 1))
        ax.set_ylabel('PC{0:d}'.format(comb[1] + 1))
        f.suptitle(area)

