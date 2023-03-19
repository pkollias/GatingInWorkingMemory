# sensory PFC
ylim = (-0.258832775674365, 1.9837228160819715)
ax = selectivity_plot(['Oscar', 'Gonzo'], 'PresentedStimulus', 'ConcatFactor',
                      ['PreDist', 'Gating', 'PostDist'], ['PreDist', 'Gating', 'PostDist'], ['PFC'],
                      ylim_arg=ylim, selective_across=True, title='Information about Presented Stimulus (Sensory)\nPFC',
                      legend_show=True, legend_cond_show=True, percent_show=False,
                      cluster_corrected=True, plot_nonzero=True, ax=None, alpha=(1, 0.3), crop=(-50, 950))
ax[0].get_figure().savefig('figures/slides/Anova_Presented_PFC_PGP.png', format='png')
# ylim = ax.get_ylim()

ax = selectivity_plot(['Oscar', 'Gonzo'], 'PresentedStimulus', 'ConcatFactor',
                      ['PreDist', 'Gating'], ['PreDist', 'Gating', 'PostDist'], ['PFC'],
                      ylim_arg=ylim, selective_across=True, title='Information about Presented Stimulus (Sensory)\nPFC',
                      legend_show=True, legend_cond_show=True, percent_show=False,
                      cluster_corrected=True, plot_nonzero=True, ax=None, alpha=(1, 0.3), crop=(-50, 950))
ax[0].get_figure().savefig('figures/slides/Anova_Presented_PFC_PG.png', format='png')

ax = selectivity_plot(['Oscar', 'Gonzo'], 'PresentedStimulus', 'ConcatFactor',
                      ['PreDist'], ['PreDist', 'Gating', 'PostDist'], ['PFC'],
                      ylim_arg=ylim, selective_across=True, title='Information about Presented Stimulus (Sensory)\nPFC',
                      legend_show=True, legend_cond_show=True, percent_show=False,
                      cluster_corrected=True, plot_nonzero=True, ax=None, alpha=(1, 0.3), crop=(-50, 950))
ax[0].get_figure().savefig('figures/slides/Anova_Presented_PFC_P.png', format='png')


# sensory IT
ylim = (-1.0546353337905356, 6.621987030858086)
ax = selectivity_plot(['Oscar', 'Gonzo'], 'PresentedStimulus', 'ConcatFactor',
                      ['PreDist', 'Gating', 'PostDist'], ['PreDist', 'Gating', 'PostDist'], ['IT'],
                      ylim_arg=ylim, selective_across=True, title='Information about Presented Stimulus (Sensory)\nIT',
                      legend_show=False, legend_cond_show=True, percent_show=False,
                      cluster_corrected=True, plot_nonzero=True, ax=None, alpha=(1, 0.3), crop=(-50, 950))
ax[0].get_figure().savefig('figures/slides/Anova_Presented_IT_PGP.png', format='png')
# ylim = ax.get_ylim()

ax = selectivity_plot(['Oscar', 'Gonzo'], 'PresentedStimulus', 'ConcatFactor',
                      ['PreDist', 'Gating'], ['PreDist', 'Gating', 'PostDist'], ['IT'],
                      ylim_arg=ylim, selective_across=True, title='Information about Presented Stimulus (Sensory)\nIT',
                      legend_show=False, legend_cond_show=True, percent_show=False,
                      cluster_corrected=True, plot_nonzero=True, ax=None, alpha=(1, 0.3), crop=(-50, 950))
ax[0].get_figure().savefig('figures/slides/Anova_Presented_IT_PG.png', format='png')

ax = selectivity_plot(['Oscar', 'Gonzo'], 'PresentedStimulus', 'ConcatFactor',
                      ['PreDist'], ['PreDist', 'Gating', 'PostDist'], ['IT'],
                      ylim_arg=ylim, selective_across=True, title='Information about Presented Stimulus (Sensory)\nIT',
                      legend_show=False, legend_cond_show=True, percent_show=False,
                      cluster_corrected=True, plot_nonzero=True, ax=None, alpha=(1, 0.3), crop=(-50, 950))
ax[0].get_figure().savefig('figures/slides/Anova_Presented_IT_P.png', format='png')


# sensory Stri
ylim = (-0.3037065344953298, 1.6670627280191153)
ax = selectivity_plot(['Oscar', 'Gonzo'], 'PresentedStimulus', 'ConcatFactor',
                      ['PreDist', 'Gating', 'PostDist'], ['PreDist', 'Gating', 'PostDist'], ['Stri'],
                      ylim_arg=ylim, selective_across=True, title='Information about Presented Stimulus (Sensory)\nStri',
                      legend_show=False, legend_cond_show=True, percent_show=False,
                      cluster_corrected=True, plot_nonzero=True, ax=None, alpha=(1, 0.3), crop=(-50, 950))
ax[0].get_figure().savefig('figures/slides/Anova_Presented_Stri_PGP.png', format='png')
# ylim = ax.get_ylim()

ax = selectivity_plot(['Oscar', 'Gonzo'], 'PresentedStimulus', 'ConcatFactor',
                      ['PreDist', 'Gating'], ['PreDist', 'Gating', 'PostDist'], ['Stri'],
                      ylim_arg=ylim, selective_across=True, title='Information about Presented Stimulus (Sensory)\nStri',
                      legend_show=False, legend_cond_show=True, percent_show=False,
                      cluster_corrected=True, plot_nonzero=True, ax=None, alpha=(1, 0.3), crop=(-50, 950))
ax[0].get_figure().savefig('figures/slides/Anova_Presented_Stri_PG.png', format='png')

ax = selectivity_plot(['Oscar', 'Gonzo'], 'PresentedStimulus', 'ConcatFactor',
                      ['PreDist'], ['PreDist', 'Gating', 'PostDist'], ['Stri'],
                      ylim_arg=ylim, selective_across=True, title='Information about Presented Stimulus (Sensory)\nStri',
                      legend_show=False, legend_cond_show=True, percent_show=False,
                      cluster_corrected=True, plot_nonzero=True, ax=None, alpha=(1, 0.3), crop=(-50, 950))
ax[0].get_figure().savefig('figures/slides/Anova_Presented_Stri_P.png', format='png')


# extended PFC
ylim = (-0.14969837984715206, 0.36930971214874475)
ax = selectivity_plot(['Oscar', 'Gonzo'], 'PresentedStimulusExtended', 'ConcatFactorExtended',
                      ['PreDist'], ['PreDist'], ['PFC'],
                      ylim_arg=ylim, selective_across=False, title='Extended Stimulus Information (PreDist)',
                      legend_show=False, legend_cond_show=True, percent_show=False,
                      cluster_corrected=True, plot_nonzero=True, ax=None, alpha=(1, 0.3), crop=(-50, 1900), ax2_y=(0, 1))
ax[0].get_figure().savefig('figures/slides/Anova_Extended_PFC_PreDist_Half.png', format='png')

ax = selectivity_plot(['Oscar', 'Gonzo'], 'PresentedStimulusExtended', 'ConcatFactorExtended',
                      ['PreDist'], ['PreDist'], ['PFC'],
                      ylim_arg=ylim, selective_across=False, title='Extended Stimulus Information (PreDist)',
                      legend_show=False, legend_cond_show=True, percent_show=False,
                      cluster_corrected=True, plot_nonzero=True, ax=None, alpha=(1, 0.3), crop=(-50, 1900), ax2_y=(0, 1))
ax = selectivity_plot(['Oscar', 'Gonzo'], 'NextStimulusExtended', 'ConcatFactorExtended',
                      ['PreDist'], ['PreDist'], ['PFC'],
                      ylim_arg=ylim, selective_across=False, title='Extended Stimulus Information (PreDist)',
                      legend_show=False, legend_cond_show=True, percent_show=False,
                      cluster_corrected=True, plot_nonzero=True, ax=ax, alpha=(0.35, 0.1), crop=(-50, 1900), ax2_y=(1, 1))
ax[0].get_figure().savefig('figures/slides/Anova_Extended_PFC_PreDist.png', format='png')

ylim = (-0.14410958260081336, 0.5219766476238107)
ax = selectivity_plot(['Oscar', 'Gonzo'], 'PresentedStimulusExtended', 'ConcatFactorExtended',
                      ['Gating'], ['Gating'], ['PFC'],
                      ylim_arg=ylim, selective_across=False, title='Extended Stimulus Information (Gating)',
                      legend_show=False, legend_cond_show=True, percent_show=False,
                      cluster_corrected=True, plot_nonzero=True, ax=None, alpha=(1, 0.3), crop=(-50, 1900), ax2_y=(0, 1))
ax = selectivity_plot(['Oscar', 'Gonzo'], 'NextStimulusExtended', 'ConcatFactorExtended',
                      ['Gating'], ['Gating'], ['PFC'],
                      ylim_arg=ylim, selective_across=False, title='Extended Stimulus Information (Gating)',
                      legend_show=False, legend_cond_show=True, percent_show=False,
                      cluster_corrected=True, plot_nonzero=True, ax=ax, alpha=(0.35, 0.1), crop=(-50, 1900), ax2_y=(1, 1))
ax[0].get_figure().savefig('figures/slides/Anova_Extended_PFC_Gating.png', format='png')

ylim = (-0.12770990132450838, 0.3507261498275955)
ax = selectivity_plot(['Oscar', 'Gonzo'], 'PresentedStimulusExtended', 'ConcatFactorExtended',
                      ['PostDist'], ['PostDist'], ['PFC'],
                      ylim_arg=ylim, selective_across=False, title='Extended Stimulus Information (PostDist)',
                      legend_show=False, legend_cond_show=True, percent_show=False,
                      cluster_corrected=True, plot_nonzero=True, ax=None, alpha=(1, 0.3), crop=(-50, 1900), ax2_y=(0, 1))
ax = selectivity_plot(['Oscar', 'Gonzo'], 'NextStimulusExtended', 'ConcatFactorExtended',
                      ['PostDist'], ['PostDist'], ['PFC'],
                      ylim_arg=ylim, selective_across=False, title='Extended Stimulus Information (PostDist)',
                      legend_show=False, legend_cond_show=True, percent_show=False,
                      cluster_corrected=True, plot_nonzero=True, ax=ax, alpha=(0.35, 0.1), crop=(-50, 1900), ax2_y=(1, 1))
ax[0].get_figure().savefig('figures/slides/Anova_Extended_PFC_PostDist.png', format='png')

# memory PFC
ylim = (-0.24516717325616877, 1.855146281742199)
ax = selectivity_plot(['Oscar', 'Gonzo'], 'GatedStimulus', 'ConcatFactor',
                      ['PreDist', 'Gating', 'PostDist', 'Target'], ['PreDist', 'Gating', 'PostDist', 'Target'], ['PFC'],
                      ylim_arg=ylim, selective_across=True, title='Information about Gated Stimulus (Memory)\nPFC',
                      legend_show=True, legend_cond_show=True, percent_show=False,
                      cluster_corrected=True, plot_nonzero=True, ax=None, alpha=(1, 0.3), crop=(-50, 950))
ax[0].get_figure().savefig('figures/slides/Anova_Gated_PFC_PGPT.png', format='png')
# ylim = ax.get_ylim()

ax = selectivity_plot(['Oscar', 'Gonzo'], 'GatedStimulus', 'ConcatFactor',
                      ['PostDist', 'Target'], ['PreDist', 'Gating', 'PostDist', 'Target'], ['PFC'],
                      ylim_arg=ylim, selective_across=True, title='Information about Gated Stimulus (Memory)\nPFC',
                      legend_show=True, legend_cond_show=True, percent_show=False,
                      cluster_corrected=True, plot_nonzero=True, ax=None, alpha=(1, 0.3), crop=(-50, 950))
ax[0].get_figure().savefig('figures/slides/Anova_Gated_PFC_PT.png', format='png')

ax = selectivity_plot(['Oscar', 'Gonzo'], 'GatedStimulus', 'ConcatFactor',
                      ['Gating', 'PostDist'], ['PreDist', 'Gating', 'PostDist', 'Target'], ['PFC'],
                      ylim_arg=ylim, selective_across=True, title='Information about Gated Stimulus (Memory)\nPFC',
                      legend_show=True, legend_cond_show=True, percent_show=False,
                      cluster_corrected=True, plot_nonzero=True, ax=None, alpha=(1, 0.3), crop=(-50, 950))
ax[0].get_figure().savefig('figures/slides/Anova_Gated_PFC_GP.png', format='png')

ax = selectivity_plot(['Oscar', 'Gonzo'], 'GatedStimulus', 'ConcatFactor',
                      ['PreDist', 'Gating'], ['PreDist', 'Gating', 'PostDist', 'Target'], ['PFC'],
                      ylim_arg=ylim, selective_across=True, title='Information about Gated Stimulus (Memory)\nPFC',
                      legend_show=True, legend_cond_show=True, percent_show=False,
                      cluster_corrected=True, plot_nonzero=True, ax=None, alpha=(1, 0.3), crop=(-50, 950))
ax[0].get_figure().savefig('figures/slides/Anova_Gated_PFC_PG.png', format='png')

ax = selectivity_plot(['Oscar', 'Gonzo'], 'GatedStimulus', 'ConcatFactor',
                      ['PreDist'], ['PreDist', 'Gating', 'PostDist', 'Target'], ['PFC'],
                      ylim_arg=ylim, selective_across=True, title='Information about Gated Stimulus (Memory)\nPFC',
                      legend_show=True, legend_cond_show=True, percent_show=False,
                      cluster_corrected=True, plot_nonzero=True, ax=None, alpha=(1, 0.3), crop=(-50, 950))
ax[0].get_figure().savefig('figures/slides/Anova_Gated_PFC_P.png', format='png')
