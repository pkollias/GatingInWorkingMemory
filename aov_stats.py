from rec_format import *
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
from itertools import groupby
import pandas as pd
import numpy as np
import warnings

def aov(df, y, x_a, x_b):

    formula = "{0} ~ C({1})".format(y, x_a) if x_b is None else '{0} ~ C({1}) + C({2}) + C({1}):C({2})'.format(y, x_a, x_b)
    model = ols(formula, df).fit()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        aov_table = anova_lm(model, typ=2)
    omega_squared(aov_table)
    return omega_squared(aov_table)


def eta_squared(aov):
    aov['eta_sq'] = 'NaN'
    aov['eta_sq'] = aov[:-1]['sum_sq']/sum(aov['sum_sq'])
    return aov


def omega_squared(aov):
    mse = aov['sum_sq'][-1]/aov['df'][-1]
    aov['omega_sq'] = 'NaN'
    aov['omega_sq'] = (aov[:-1]['sum_sq']-(aov[:-1]['df']*mse))/(sum(aov['sum_sq'])+mse)
    return aov

# TODO: when saving omega_sq value round up to x-th digit

def aov_2_shuffles(df, y, x_a, x_b, x_ab, num_shuffles=1, group_column_list=[]):

    def shuffle_df(df, y, x_a, x_b, shuffle_i):
        return pd.concat([df[x_a].sample(frac=1, random_state=shuffle_i).reset_index(drop=True),
                          df[x_b].sample(frac=1, random_state=shuffle_i).reset_index(drop=True),
                          df[y].reset_index(drop=True)], axis=1)

    observed_aov = aov(df, y, x_a, x_b)

    # anovas
    anovas = pd.DataFrame(columns=[x_a, x_b, x_ab])
    # observed
    anovas.loc['observed'] = list(observed_aov['omega_sq'][0:3])
    # shuffled
    for shuffle_i in range(num_shuffles - 1):
        if group_column_list == []:
            df_shuffle = shuffle_df(df, y, x_a, x_b, shuffle_i)
        else:
            df_grouper = df.groupby(group_column_list)
            dfg_list = [df_grouper.get_group(k) for k in df_grouper.groups.keys()]
            df_shuffle_list = [shuffle_df(dfg, y, x_a, x_b, shuffle_i)
                               for dfg in dfg_list]
            df_shuffle = pd.concat(df_shuffle_list, axis=0, ignore_index=True)
        anovas.loc['shuffle_{0:04d}'.format(shuffle_i)] = list(aov(df_shuffle, y, x_a, x_b)['omega_sq'][0:3])

    return (anovas, observed_aov)


def aov_1_shuffles(df, y, x, num_shuffles=1, group_column_list=[]):

    def shuffle_df(df, y, x, shuffle_i):
        return pd.concat([df[x].sample(frac=1, random_state=shuffle_i).reset_index(drop=True),
                          df[y].reset_index(drop=True)], axis=1)

    observed_aov = aov(df, y, x, None)

    # anovas
    anovas = pd.DataFrame(columns=[x])
    # observed
    anovas.loc['observed'] = [observed_aov['omega_sq'][0]]
    # shuffled
    for shuffle_i in range(num_shuffles - 1):
        if group_column_list == []:
            df_shuffle = shuffle_df(df, y, x, shuffle_i)
        else:
            df_grouper = df.groupby(group_column_list)
            dfg_list = [df_grouper.get_group(k) for k in df_grouper.groups.keys()]
            df_shuffle_list = [shuffle_df(dfg, y, x, shuffle_i)
                               for dfg in dfg_list]
            df_shuffle = pd.concat(df_shuffle_list, axis=0, ignore_index=True)
        anovas.loc['shuffle_{0:04d}'.format(shuffle_i)] = [aov(df_shuffle, y, x, None)['omega_sq'][0]]

    return (anovas, observed_aov)


def aov_2_shuffle_results(df, anovas, observed_aov, y, x_a, x_b, x_ab):

    unit_results = {'means': {x_a: df[[x_a, y]].groupby(x_a).mean().dropna(),
                              x_b: df[[x_b, y]].groupby(x_b).mean().dropna(),
                              x_ab.format(x_a, x_b): df[[x_a, x_b, y]].groupby([x_a, x_b]).mean().dropna()},
                    'anova': observed_aov,
                    'shuffles': anovas}

    return unit_results


def aov_1_shuffle_results(df, anovas, observed_aov, y, x):

    unit_results = {'means': {x: df[[x, y]].groupby(x).mean().dropna()},
                    'anova': observed_aov,
                    'shuffles': anovas}

    return unit_results


def aov_shuffle_and_results(df, y, x_list, num_shuffles, group_column_list = []):

    if len(x_list) == 2:
        x_a = x_list[0]
        x_b = x_list[1]
        x_ab = interaction_term(x_a, x_b)
        anovas, observed_aov = aov_2_shuffles(df, y, x_a, x_b, x_ab, num_shuffles, group_column_list)
        unit_results = aov_2_shuffle_results(df, anovas, observed_aov, y, x_a, x_b, x_ab)
    elif len(x_list) == 1:
        x = x_list[0]
        anovas, observed_aov = aov_1_shuffles(df, y, x, num_shuffles, group_column_list)
        unit_results = aov_1_shuffle_results(df, anovas, observed_aov, y, x)

    return unit_results


def sorted_percentile(sorted_array, perc):

    pos = (len(sorted_array) - 1) * perc / 100
    x1 = int(np.floor(pos))
    x2 = x1 + 1 if x1 == pos else int(np.ceil(pos))
    y1 = sorted_array[x1]
    y2 = sorted_array[x2]
    slope = (y2 - y1) / (x2 - x1)
    intercept = y1 - slope * x1
    return slope * pos + intercept

# round omega squared values and correct for zero shuffles (all values equal to zero)
# return unit validity, omega distribution,  threshold, omega difference above threshold, and zscored observed
def evaluate_anova_shuffles(shuffles_series, percentile, shuffles):

    rounding_decimal = 10
    anova_round = lambda x: round(x, rounding_decimal)
    shuffles_index = [shuffle_to_name(ii) for ii in range(shuffles)]

    if len(shuffles_series) == 1:
        valid = False
        omega_sq_distr = pd.Series(data=0, index=shuffles_index)
        omega_sq_observed = omega_sq_distr.loc[shuffle_to_name(0)]
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


def clusters_from_shuffle_list(omega_sq_diff_list):
    # returns cluster list of cluster_result cr where each cluster consists of a list of consecutive significant timebin
    # tuples of the time value and the cluster difference
    clusters_list = [list(gr) for val, gr in groupby(omega_sq_diff_list, key=lambda x: x[1] > 0) if val]
    # convert that into a list of tuples [(timebins, diff_values), ...]
    return [tuple(map(list, zip(*cl))) for cl in clusters_list]


def get_clusters(clusters):

    return [{'bins': cv[0], 'val': sum(cv[1])} for cv in clusters]
