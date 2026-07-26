###################################################################################################
# SERVER
###################################################################################################

server <- function(input, output, session) {
  
  page_overview_server(input, output, session, annotated_genome, gem_network)
  
  page_genome_server(input, output, session, annotated_genome)
  page_phac_server(input, output, session)
  page_gem_server(input, output, session, gem_network)
  page_fba_server(input, output, session, fba_results, pfba_flux_distribution, production_envelope)
  page_flexibility_server(input, output, session, pathway_flexibility)
  
  page_environmental_screening_server(
    input, output, session,
    medium_screen_results, minimal_medium_exchanges, aeration_carbon_results
  )
  
  page_robustness_server(
    input, output, session,
    robust_1d_glucose, robust_1d_maltose, robust_1d_nitrogen, robust_1d_oxygen,
    robust_2d_glc_nh4, robust_2d_malt_nh4, robust_2d_glc_o2
  )
  
  page_phpp_server(
    input, output, session,
    phpp_glucose_grid, phpp_glucose_transitions,
    phpp_maltose_grid, phpp_maltose_transitions
  )
  
  page_dfba_server(
    input, output, session,
    dfba_glucose_standard, dfba_glucose_nlimited, dfba_glucose_twostage,
    dfba_maltose_standard, dfba_maltose_nlimited, dfba_maltose_twostage
  )
  
  page_diagnostics_server(
    input, output, session,
    diag_shadow_price_growth, diag_reduced_cost_growth,
    diag_shadow_price_phb, diag_reduced_cost_phb
  )
  
  page_bottlenecks_server(input, output, session, bottleneck_results)
  
  page_multifva_server(
    input, output, session,
    multifva_100_standard, multifva_100_loopless,
    multifva_90_standard,  multifva_90_loopless,
    multifva_50_standard,  multifva_50_loopless,
    pathway_flexibility
  )
  
  page_strain_design_server(
    input, output, session,
    knockout_essentiality, knockout_single, knockout_double
  )
  
  page_cost_server(input, output, session)
  
}
