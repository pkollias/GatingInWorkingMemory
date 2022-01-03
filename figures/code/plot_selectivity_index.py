from rec import *
from aov_stats import *
import matplotlib.pyplot as plt
from sklearn.linear_model import (LinearRegression, HuberRegressor)
import seaborn as sns
from scipy.odr import *
from itertools import chain
from matplotlib.gridspec import GridSpec
from versioning import *

def main():

    version_aov = 'GeneralizedGating'

    # init tables and slices
    md = MetaData()
    db = md.db_base_loader(['units', 'physiology'])
    units, physiology = db['units'], db['physiology']


    ####################################################

    def anova_significance(physiology_dict_entry):
        return bool(physiology_dict_entry['GatingCondSpecialized']['significant'])

    def point_from_line_distance(p1, p2, p3):
        np1, np2, np3 = np.array(p1), np.array(p2), np.array(p3)
        return np.linalg.norm(np.cross(np2-np1, np1-np3))/np.linalg.norm(np2-np1)

    point_from_diagonal_distance = lambda p: point_from_line_distance((0, 0), (1, 1), p) * np.sign(p[1] - p[0])

    def r_squared(x, y, model):
        yhat = model.predict(x)
        SS_Residual = sum((y - yhat) ** 2)
        SS_Total = sum((y - np.mean(y)) ** 2)
        return 1 - (float(SS_Residual)) / SS_Total

    def predict(x, coef, intercept):
        return coef * np.array(x) + intercept

    def linear_func(p, x):
        m, c = p
        return m * x + c

    def plot_modulation_scatter(fig, ax, x, y, coef, intercept, reg, coef_inv, intercept_inv, reg_inv, coef_o, intercept_o, mode, fits):
        ax.scatter(x, y, c=color[area], s=30, linewidths=linewidths, edgecolors='black')
        ax.axis('square')
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        minlim = min(xlim[0], ylim[0])
        maxlim = min(xlim[1], ylim[1])
        ax.text(.9, .95, 'y = x', horizontalalignment='right', verticalalignment='top', transform=ax.transAxes, color='black')

        # if fits:
        #     ax.plot(xlim, reg.predict(np.array(xlim).reshape(-1, 1)), color='red')
        #     ax.plot(xlim, predict(np.array(xlim).reshape(-1, 1), coef_inv, intercept_inv), color='green')
        #     ax.plot(xlim, predict(np.array(xlim).reshape(-1, 1), coef_o, intercept_o), color='blue')
        #     regression_str = 'y = {0:.2f}*x {1:+.2f}\n$R^2 = {2:.2f}$ (y~x)'.format(coef, intercept, reg.score(x, y))
        #     ax.text(.95, 0.05, regression_str, horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes, color='red')
        #     regression_inv_str = 'y = {0:.2f}*x {1:+.2f}\n$R^2 = {2:.2f}$ (x~y)'.format(coef_inv, intercept_inv, reg_inv.score(y, x))
        #     ax.text(.95, 0.25, regression_inv_str, horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes, color='green')
        #     regression_o_str = 'y = {0:.2f}*x {1:+.2f} (xy diag)'.format(coef_o, intercept_o)
        #     ax.text(.95, 0.45, regression_o_str, horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes, color='blue')

        ax.plot([minlim, maxlim], [minlim, maxlim], color='k')
        center = 1 if mode == 'Modulation' else 0
        ax.plot(list(xlim), [center, center], color='k', ls=':')
        ax.plot([center, center], list(ylim), color='k', ls=':')
        return (fig, ax)

    ####################################################

    # Gating modulation scatter plot
    color_list = ['#d817a5', '#48b213', '#727272']
    line_list = ['PFC', 'Stri', 'IT']
    color = dict(zip(line_list, color_list))

    for version_fr in ['WindowSampleShift', 'WindowDelayShift']:

        v_fr_params = version_fr_params(version_fr)
        t_start = v_fr_params['t_start']
        t_end = v_fr_params['t_end']
        timebin = v_fr_params['timebin']
        timestep = v_fr_params['timestep']
        src_filename = md.proc_dest_path(path.join('BehavioralUnits', 'Anova', version_aov,
                                                      behunit_params_str(version_fr, timebin, timestep, t_start, t_end),
                                                      'selectivity_index'),
                                            'physiology_df.pkl')
        physiology_df = md.np_loader(src_filename)
        df = pd.concat([units, physiology_df], axis=1)
        df = df.loc[:, ~df.columns.duplicated()]
        df = df.loc[df['UnitNum'].ne(0) & df['RatingCode'].ne(7)]

        src_filename = md.proc_dest_path(path.join('BehavioralUnits', 'Anova', version_aov,
                                                   behunit_params_str(version_fr, timebin, timestep, t_start, t_end),
                                                   'summarize'),
                                         'physiology_dict.pkl')
        physiology_dict = md.np_loader(src_filename)

        for x_label, y_label,\
            mode, stim_group, scale in [(['PreDistModulation_All'], ['GatingModulation_All'], 'Modulation', 'Stim_Average', 'normal'),
                                        (['PreDistModulation_All'], ['GatingModulation_All'], 'Modulation', 'Stim_Average', 'log'),
                                        (['PreDistDifference_All'], ['GatingDifference_All'], 'Difference', 'Stim_Average', 'normal'),
                                        (['PreDistContrast_All'], ['GatingContrast_All'], 'Contrast', 'Stim_Average', 'normal'),
                                        (['PreDistModulation_S11', 'PreDistModulation_S12', 'PreDistModulation_S21', 'PreDistModulation_S22'],
                                         ['GatingModulation_S11', 'GatingModulation_S12', 'GatingModulation_S21', 'GatingModulation_S22'], 'Modulation', 'Stim_Individual', 'normal'),
                                        (['PreDistModulation_S11', 'PreDistModulation_S12', 'PreDistModulation_S21', 'PreDistModulation_S22'],
                                         ['GatingModulation_S11', 'GatingModulation_S12', 'GatingModulation_S21', 'GatingModulation_S22'], 'Modulation', 'Stim_Individual', 'log'),
                                        (['PreDistDifference_S11', 'PreDistDifference_S12', 'PreDistDifference_S21', 'PreDistDifference_S22'],
                                         ['GatingDifference_S11', 'GatingDifference_S12', 'GatingDifference_S21', 'GatingDifference_S22'], 'Difference', 'Stim_Individual', 'normal'),
                                        (['PreDistContrast_S11', 'PreDistContrast_S12', 'PreDistContrast_S21', 'PreDistContrast_S22'],
                                         ['GatingContrast_S11', 'GatingContrast_S12', 'GatingContrast_S21', 'GatingContrast_S22'], 'Contrast', 'Stim_Individual', 'normal')]:
            for regressor in [LinearRegression, HuberRegressor]:

                fig = plt.figure(constrained_layout=True, figsize=(16, 9))
                gs = GridSpec(6, 9, figure=fig)
                fig.suptitle('{1:s}\nFiring Rate {2:s}\n{3:s}\n{4:s}\n{5:s}-scale axes'.format(_, version_fr, mode, regressor.__name__, stim_group, scale))

                for area_index, area in enumerate(['PFC', 'Stri', 'IT']):


                    df_area = df.loc[df['Area'].eq(area)]
                    df_area_filtered = df_area.loc[df_area[x_label + y_label].lt(100).all(axis=1)].dropna()
                    if scale == 'log':
                        df_area_filtered = df_area_filtered.drop(df_area_filtered.loc[df_area_filtered[x_label + y_label].eq(0).any(axis=1)].index)
                    significance_mask = list(chain.from_iterable([[anova_significance(physiology_dict['All'][key])
                                                                   for key in list(df_area_filtered.index)]
                                                                  for _ in range(len(x_label))]))
                    linewidths = [int(significance) * 2 for significance in significance_mask]
                    sign_index = [np.where(significance_mask)[0], np.where(np.invert(significance_mask))[0]]


                    x = np.array(df_area_filtered[x_label]).reshape(-1, 1, order='F')
                    y = np.array(df_area_filtered[y_label]).reshape(-1, 1, order='F')
                    # fit inverse relationship and then inverse

                    fit_x = np.log(x) if scale == 'log' else x
                    fit_y = np.log(y) if scale == 'log' else y
                    if regressor == LinearRegression:
                        reg = regressor()
                        reg.fit(fit_x, fit_y)
                        coef = reg.coef_[0][0]
                        intercept = reg.intercept_[0]
                    else:
                        reg = regressor(epsilon=2.5)
                        reg.fit(fit_x, fit_y)
                        coef = reg.coef_[0]
                        intercept = reg.intercept_

                    if regressor == LinearRegression:
                        reg_inv = regressor()
                        reg_inv.fit(fit_y, fit_x)
                        coef_inv = 1/reg_inv.coef_[0][0]
                        intercept_inv = -reg_inv.intercept_[0] * coef_inv
                    else:
                        reg_inv = regressor(epsilon=2.5)
                        reg_inv.fit(fit_y, fit_x)
                        coef_inv = 1/reg_inv.coef_[0]
                        intercept_inv = -reg.intercept_ * coef_inv


                    if regressor == HuberRegressor:
                        non_outliers_mask = np.invert(reg.outliers_)
                        x = x[non_outliers_mask]
                        y = y[non_outliers_mask]
                        significance_mask = list(np.array(significance_mask)[non_outliers_mask])
                        linewidths = list(np.array(linewidths)[non_outliers_mask])



                    # create dataframe
                    coords_list = list(zip(list(x.reshape(-1)), list(y.reshape(-1))))
                    residuals_from_diagonal = [point_from_diagonal_distance(p) for p in coords_list]
                    residuals_df = pd.DataFrame({'Residuals': residuals_from_diagonal, 'Significant': significance_mask})


                    # Diagonal Squares Regression
                    # Create a model for fitting.
                    linear_model = Model(linear_func)
                    # Create a RealData object using our initiated data from above.
                    data = RealData(fit_x.transpose()[0], fit_y.transpose()[0])
                    # Set up ODR with the model and data.
                    odr = ODR(data, linear_model, beta0=[1., 0.])
                    # Run the regression.
                    out = odr.run()
                    # Use the in-built pprint method to give us results.
                    coef_o = out.beta[0]
                    intercept_o = out.beta[1]

                    ax1 = fig.add_subplot(gs[:3, slice(3 * area_index, 3 * (area_index + 1))])
                    ax1a = fig.add_subplot(gs[3, slice(3 * area_index, 3 * (area_index + 1))])
                    ax1b = fig.add_subplot(gs[4, slice(3 * area_index, 3 * (area_index + 1))])
                    ax1c = fig.add_subplot(gs[5, slice(3 * area_index, 3 * (area_index + 1))])

                    params = (fig, ax1, x, y, coef, intercept, reg, coef_inv, intercept_inv, reg_inv, coef_o, intercept_o, mode)
                    if scale == 'log':
                        ax1.set_xscale('log')
                        ax1.set_yscale('log')
                        fig, ax = plot_modulation_scatter(*params, True)
                    else:
                        fig, ax = plot_modulation_scatter(*params, True)
                    mod_char = '/' if mode == 'Modulation' else '-'
                    numer_x_str = '($FR_{predist}$' + mod_char + '$FR_{fixation}$)'
                    denom_x_str = ' / ($FR_{predist} + FR_{fixation}$)' if mode == 'Contrast' else ''
                    numer_y_str = ' ($FR_{gating}$' + mod_char + '$FR_{fixation}$)'
                    denom_y_str = ' / ($FR_{gating} + FR_{fixation}$)' if mode == 'Contrast' else ''
                    ax1.set_xlabel(numer_x_str + denom_x_str)
                    ax1.set_ylabel(numer_y_str + denom_y_str)
                    ax1.set_title('{0:s}'.format(area))

                    sns.histplot(data=residuals_df, x="Residuals", hue="Significant", element='step', fill=True, stat='probability', common_norm=False, ax=ax1a, alpha=0.075)
                    ax1a.vlines(0, ax1a.get_ylim()[0], ax1a.get_ylim()[1], color='k')
                    ax1a.set_xlabel('Distance')
                    ax1a.set_ylabel('counts')
                    ax1a.legend(['Significant $\omega^2$', 'Non-significant'])
                    ax1a.set_title('Distance from Diagonal Histogram')

                    sns.histplot(data=pd.DataFrame({'x': x.flatten(), 'Significant': significance_mask}), x="x", hue="Significant", element='step', fill=True, stat='probability', common_norm=False, ax=ax1b, alpha=0.075)
                    ax1b.vlines(0, ax1b.get_ylim()[0], ax1b.get_ylim()[1], color='k')
                    ax1b.set_xlim(ax1.get_xlim())
                    ax1b.set_xlabel('PreDist')
                    ax1b.set_ylabel('counts')
                    ax1b.legend(['Significant $\omega^2$', 'Non-significant'])

                    sns.histplot(data=pd.DataFrame({'y': y.flatten(), 'Significant': significance_mask}), x="y", hue="Significant", element='step', fill=True, stat='probability', common_norm=False, ax=ax1c, alpha=0.075)
                    ax1c.vlines(0, ax1c.get_ylim()[0], ax1c.get_ylim()[1], color='k')
                    ax1c.set_xlim(ax1.get_ylim())
                    ax1c.set_xlabel('Gating')
                    ax1c.set_ylabel('counts')
                    ax1c.legend(['Significant $\omega^2$', 'Non-significant'])


                plt.savefig('Modulation_{1:s}_{2:s}_{3:s}_{4:s}_{5:s}.png'.format(_, version_fr, mode, regressor.__name__, stim_group, scale), format='png')
                plt.close()


        # Depth of selectivity cum hist
        for period in ['Overall', 'Distractor', 'Gating']:

            area_list = ['PFC', 'Stri', 'IT']

            fig = plt.figure()
            ax2 = plt.gca()
            period_str = 'DepthOfSelectivityIndex_{0:s}'.format(period)
            df_drop = df.drop(df.loc[df[period_str].gt(1)].index)[['Area', period_str]].dropna()
            for area in area_list:
                df_area = df_drop.loc[df_drop['Area'].eq(area)]
                sns.histplot(df_area, x=period_str, color=color[area], stat='probability', element='poly', bins=30, cumulative=True, fill=False)
            ax2.set_xlabel('Depth of Stimulus Selectivity Index')
            ax2.set_ylabel('Cumulative Distribution of Population')
            ax2.legend(area_list)
            ax2.set_title('{0:s} {1:s}'.format(period, version_fr))
            plt.savefig('DoSSI_{0:s}_{1:s}_Cumul.png'.format(period, version_fr), format='png')
            plt.close()

            fig = plt.figure()
            ax3 = plt.gca()
            period_str = 'DepthOfSelectivityIndex_{0:s}'.format(period)
            df_drop = df.drop(df.loc[df[period_str].gt(1)].index)[['Area', period_str]].dropna()
            for area in area_list:
                df_area = df_drop.loc[df_drop['Area'].eq(area)]
                sns.histplot(df_area, x=period_str, color=color[area], stat='probability', element='poly', bins=30, cumulative=False, fill=False)
            ax3.set_xlabel('Depth of Stimulus Selectivity Index')
            ax3.set_ylabel('Probability Distribution of Population')
            ax3.legend(area_list)
            ax3.set_title('{0:s} {1:s}'.format(period, version_fr))
            plt.savefig('DoSSI_{0:s}_{1:s}_Prob.png'.format(period, version_fr), format='png')
            plt.close()

            fig = plt.figure()
            ax4 = plt.gca()
            period_str = 'StimulusSelectivityIndex_{0:s}'.format(period)
            df_drop = df.loc[df[period_str].lt(8)].dropna()
            for area in area_list:
                df_area = df_drop.loc[df_drop['Area'].eq(area)]
                sns.histplot(df_area, x=period_str, color=color[area], stat='probability', element='poly', bins=30, cumulative=True, fill=False)
            ax4.set_xlabel('Stimulus Selectivity Index')
            ax4.set_ylabel('Cumulative Distribution of Population')
            ax4.legend(area_list)
            ax4.set_title('{0:s} {1:s}'.format(period, version_fr))
            plt.savefig('SSI_{0:s}_{1:s}_Cumul.png'.format(period, version_fr), format='png')
            plt.close()

            fig = plt.figure()
            ax5 = plt.gca()
            period_str = 'StimulusSelectivityIndex_{0:s}'.format(period)
            df_drop = df.loc[df[period_str].lt(8)].dropna()
            for area in area_list:
                df_area = df_drop.loc[df_drop['Area'].eq(area)]
                sns.histplot(df_area, x=period_str, color=color[area], stat='probability', element='poly', bins=30, cumulative=False, fill=False)
            ax5.set_xlabel('Stimulus Selectivity Index')
            ax5.set_ylabel('Probability Distribution of Population')
            ax5.legend(area_list)
            ax5.set_title('{0:s} {1:s}'.format(period, version_fr))
            plt.savefig('SSI_{0:s}_{1:s}_Prob.png'.format(period, version_fr), format='png')
            plt.close()




main()
