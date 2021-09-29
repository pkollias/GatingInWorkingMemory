
# plot script
ax = selectivity_plot(['Gonzo', 'Oscar'], 'PresentedStimulus', 'ConcatFactor', ['PreDist', 'Gating', 'PostDist'], ['PreDist', 'Gating', 'PostDist'], ['PFC'], (), None, False, '', True, True)
ax.get_figure().savefig('figures/slides/PresentedPFC.png', format='png')

ax = selectivity_plot(['Gonzo', 'Oscar'], 'PresentedStimulus', 'ConcatFactor', ['PreDist', 'Gating', 'PostDist'], ['PreDist', 'Gating', 'PostDist'], ['Stri'], (), None, False, '', True, False)
ax.get_figure().savefig('figures/slides/PresentedStri.png', format='png')

ax = selectivity_plot(['Gonzo', 'Oscar'], 'PresentedStimulus', 'ConcatFactor', ['PreDist', 'Gating', 'PostDist'], ['PreDist', 'Gating', 'PostDist'], ['IT'], (), None, False, '', True, False)
ax.get_figure().savefig('figures/slides/PresentedIT.png', format='png')

ax = selectivity_plot(['Gonzo', 'Oscar'], 'GatedStimulus', 'ConcatFactor', ['Gating', 'PostDist', 'Target'], ['Gating', 'PostDist', 'Target'], ['PFC'], (), None, False, '', True, True)
ax.get_figure().savefig('figures/slides/GatedPFC.png', format='png')

ax = dpca_pev_plot('PostStimulusRuleStim', 'ConcatFactor2', 'Balance', 'PFC', 'Gonzo_Oscar', 20, 3, inner_show=False, legend_show=True, title_str='PFC')
ax.get_figure().savefig('figures/slides/SensMemPEVPFC_outer_nolabels.png', format='png')

ax = dpca_pev_plot('PostStimulusRuleStim', 'ConcatFactor2', 'Balance', 'PFC', 'Gonzo_Oscar', 20, 3, inner_show=False, legend_show=True, title_str='PFC')
ax.get_figure().savefig('figures/slides/SensMemPEVPFC_outer.png', format='png')

ax1, ax2 = dpca_state_space_plot('PostStimulusRuleStim', 'ConcatFactor2', 'Balance', 'PFC', 'Gonzo_Oscar', 20, [('x', -1), ('t', 0)], condition_list_slice=False, legend_show=True, azim=-60, elev=30, legend_letters=('s', 'm'))
ax1.get_figure().savefig('figures/slides/SensMem_PFC_xtPC1.png', format='png')

ax1, ax2 = dpca_state_space_plot('PostStimulusRuleStim', 'ConcatFactor2', 'Balance', 'PFC', 'Gonzo_Oscar', 20, [('t', 0), ('t', 1)], condition_list_slice=False, legend_show=False, azim=-60, elev=30, legend_letters=('s', 'm'))
# get ani

# PFC - SensMem
for version_factor in ['GatPostStimulusRuleStim', 'PostStimulusRuleStim', 'StimulugGating', 'StimulusGatingBool', 'RuleStimGatingNull']:
    for area in ['PFC', 'Stri', 'IT']:
        for subject in ['Gonzo_Oscar', 'Gonzo', 'Oscar']:
            for m in ['t', 's', 'm', 'g', 'sm', 'sg', 'mg']:
                ind = 0
                for margin_dim_list in [[('x', -1), (m, 0)], [('x', -1), (m, 1)], [('x', -1), (m, 2)],
                                        [(m, 0), (m, 1)], [(m, 0), (m, 2)], [(m, 0), (m, 1), (m, 2)]]:
                    try:
                        ind = ind + 1
                        ax1, ax2 = dpca_state_space_plot(version_factor, 'ConcatFactor2', 'Balance', area, 'Gonzo_Oscar', 20,
                                                     margin_dim_list, condition_list_slice=False, legend_show=False, azim=123,
                                                     elev=30, legend_letters=('s', 'g'))
                    except:
                        pass



ax1, ax2 = dpca_state_space_plot('PostStimulusRuleStim', 'ConcatFactor2', 'Balance', 'PFC', 'Gonzo_Oscar', 20, [('sm', 0), ('sm', 1)], condition_list_slice=False, legend_show=False, azim=-9.74, elev=49.63, legend_letters=('s', 'm'))
# get ani

ax1, ax2 = dpca_state_space_plot('PostStimulusRuleStim', 'ConcatFactor2', 'Balance', 'Stri', 'Gonzo_Oscar', 20, [('sm', 0), ('sm', 1)], condition_list_slice=False, legend_show=False, azim=-9.74, elev=49.63, legend_letters=('s', 'm'))
# get ani

ax1, ax2 = dpca_state_space_plot('PostStimulusRuleStim', 'ConcatFactor2', 'Balance', 'IT', 'Gonzo_Oscar', 20, [('s', 0), ('s', 1)], condition_list_slice=False, legend_show=True, azim=-9.74, elev=49.63, legend_letters=('s', 'm'))
ax1.get_figure().savefig('figures/slides/SensMem_IT_sPC1sPC2.png', format='png')

ax1, ax2 = dpca_state_space_plot('PostStimulusRuleStim', 'ConcatFactor2', 'Balance', 'IT', 'Gonzo_Oscar', 20, [('s', 0), ('s', 2)], condition_list_slice=False, legend_show=False, azim=-9.74, elev=49.63, legend_letters=('s', 'm'))
ax1.get_figure().savefig('figures/slides/SensMem_IT_sPC1sPC3.png', format='png')

ax1, ax2 = dpca_state_space_plot('PostStimulusRuleStim', 'ConcatFactor2', 'Balance', 'IT', 'Gonzo_Oscar', 20, [('x', -1), ('t', 0)], condition_list_slice=False, legend_show=False, azim=-9.74, elev=49.63, legend_letters=('s', 'm'))
ax1.get_figure().savefig('figures/slides/SensMem_IT_xtPC1.png', format='png')

ax1, ax2 = dpca_state_space_plot('PostStimulusRuleStim', 'ConcatFactor2', 'Balance', 'PFC', 'Gonzo_Oscar', 20, [('m', 0), ('m', 1)], condition_list_slice=False, legend_show=True, azim=-9.74, elev=49.63, legend_letters=('s', 'm'))
ax1.get_figure().savefig('figures/slides/SensMem_PFC_mPC1mPC2.png', format='png')

ax1, ax2 = dpca_state_space_plot('PostStimulusRuleStim', 'ConcatFactor2', 'Balance', 'IT', 'Gonzo_Oscar', 20, [('m', 0), ('m', 1)], condition_list_slice=False, legend_show=False, azim=-9.74, elev=49.63, legend_letters=('s', 'm'))
ax1.get_figure().savefig('figures/slides/SensMem_IT_mPC1mPC2.png', format='png')

ax1, ax2 = dpca_state_space_plot('PostStimulusRuleStim', 'ConcatFactor2', 'Balance', 'Stri', 'Gonzo_Oscar', 20, [('m', 0), ('m', 1)], condition_list_slice=False, legend_show=False, azim=-9.74, elev=49.63, legend_letters=('s', 'm'))
ax1.get_figure().savefig('figures/slides/SensMem_Stri_mPC1mPC2.png', format='png')


ax = dpca_pev_plot('StimulusGating', 'ConcatFactor2', 'Balance', 'PFC', 'Gonzo_Oscar', 20, 3, inner_show=True, legend_show=True, title_str='PFC')
ax.get_figure().savefig('figures/slides/SensGatPEVPFC.png', format='png')

ax = dpca_pev_plot('StimulusGating', 'ConcatFactor2', 'Balance', 'IT', 'Gonzo_Oscar', 20, 3, inner_show=True, legend_show=False, title_str='IT')
ax.get_figure().savefig('figures/slides/SensGatPEVIT.png', format='png')

ax = dpca_pev_plot('StimulusGating', 'ConcatFactor2', 'Balance', 'Stri', 'Gonzo_Oscar', 20, 3, inner_show=True, legend_show=False, title_str='Stri')
ax.get_figure().savefig('figures/slides/SensGatPEVStri.png', format='png')

ax1, ax2 = dpca_state_space_plot('StimulusGating', 'ConcatFactor2', 'Balance', 'PFC', 'Gonzo_Oscar', 20, [('g', 0), ('g', 1), ('g', 2)], condition_list_slice=False, legend_show=True, azim=-132, elev=36, legend_letters=('', ''))
ax1.get_figure().savefig('figures/slides/SensGat_PFC_gPC1gPC2gPC3.png', format='png')

ax1, ax2 = dpca_state_space_plot('StimulusGating', 'ConcatFactor2', 'Balance', 'PFC', 'Gonzo_Oscar', 20, [('g', 0), ('g', 1), ('g', 2)], condition_list_slice=False, legend_show=False, azim=-139, elev=22, legend_letters=('', ''))
# grab ani

ax1, ax2 = dpca_state_space_plot('StimulusGating', 'ConcatFactor2', 'Balance', 'Stri', 'Gonzo_Oscar', 20, [('g', 0), ('g', 1), ('g', 2)], condition_list_slice=False, legend_show=False, azim=-118, elev=33, legend_letters=('', ''))
ax1.get_figure().savefig('figures/slides/SensGat_Stri_gPC1gPC2gPC3.png', format='png')

ax = dpca_pev_plot('StimulusGatingBool', 'ConcatFactor2', 'Balance', 'PFC', 'Gonzo_Oscar', 20, 3, inner_show=True, legend_show=True, title_str='PFC')
ax.get_figure().savefig('figures/slides/SensGatBoolPEVPFC.png', format='png')

ax = dpca_pev_plot('StimulusGatingBool', 'ConcatFactor2', 'Balance', 'IT', 'Gonzo_Oscar', 20, 3, inner_show=True, legend_show=False, title_str='IT')
ax.get_figure().savefig('figures/slides/SensGatBoolPEVIT.png', format='png')

ax = dpca_pev_plot('StimulusGatingBool', 'ConcatFactor2', 'Balance', 'Stri', 'Gonzo_Oscar', 20, 3, inner_show=True, legend_show=False, title_str='Stri')
ax.get_figure().savefig('figures/slides/SensGatBoolPEVStri.png', format='png')

ax1, ax2 = dpca_state_space_plot('StimulusGatingBool', 'ConcatFactor2', 'Balance', 'Stri', 'Gonzo_Oscar', 20, [('x', -1), ('g', 0)], condition_list_slice=False, legend_show=True, azim=-118, elev=33, legend_letters=('', ''))
ax1.get_figure().savefig('figures/slides/SensGatBool_Stri_xtgPC1.png', format='png')