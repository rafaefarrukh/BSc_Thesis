###################################################################################################
# SERVER: Shadow Prices & Reduced Costs
###################################################################################################

page_diagnostics_server <- function(input, output, session,
                                     diag_shadow_price_growth, diag_reduced_cost_growth,
                                     diag_shadow_price_phb,    diag_reduced_cost_phb) {

  # Joins a metabolite-level shadow-price table with its corresponding exchange-reaction-level
  # reduced-cost table. Every exchanged metabolite `X_e` has a matching exchange reaction
  # `EX_X_e`, so the join key is derived by stripping the "EX_" prefix from the reaction id.
  build_diag_table <- function(shadow_df, reduced_df) {
    reduced_join <- reduced_df %>%
      mutate(metabolite_id = sub("^EX_", "", reaction_id)) %>%
      select(metabolite_id, reduced_cost, lower_bound, upper_bound)

    shadow_df %>%
      left_join(reduced_join, by = "metabolite_id") %>%
      arrange(shadow_price) %>%
      transmute(
        `Metabolite`            = metabolite_name,
        `Metabolite ID`         = metabolite_id,
        `Shadow Price`          = signif(shadow_price, 4),
        `Exchange Reduced Cost` = signif(reduced_cost, 4),
        `Lower Bound`           = lower_bound,
        `Upper Bound`           = upper_bound
      )
  }

  output$diag_growth_opt_table <- renderDT({
    req(diag_shadow_price_growth, diag_reduced_cost_growth)

    build_diag_table(diag_shadow_price_growth, diag_reduced_cost_growth) %>%
      datatable(
        rownames  = FALSE,
        selection = "single",
        options   = list(pageLength = 15, scrollX = TRUE)
      )
  })

  output$diag_phb_opt_table <- renderDT({
    req(diag_shadow_price_phb, diag_reduced_cost_phb)

    build_diag_table(diag_shadow_price_phb, diag_reduced_cost_phb) %>%
      datatable(
        rownames  = FALSE,
        selection = "single",
        options   = list(pageLength = 15, scrollX = TRUE)
      )
  })

}

###################################################################################################
