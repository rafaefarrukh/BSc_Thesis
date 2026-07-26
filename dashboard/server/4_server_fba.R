###################################################################################################
# SERVER: FBA & Production Envelope
###################################################################################################

page_fba_server <- function(input, output, session, fba_results, pfba_flux_distribution, production_envelope) {
  
  # FBA SCENARIOS TABLE ############################
  output$fba_scenarios_table <- renderDT({
    req(fba_results)
    
    fba_results %>%
      transmute(
        `Scenario` = scenario,
        `Status` = status,
        `Growth Rate (1/h)` = round(growth_rate_1_h, 4),
        `PHB Flux (mmol/gDW/h)` = round(phb_flux, 4)
      ) %>%
      datatable(
        rownames = FALSE,
        selection = "single",
        options = list(
          paging = FALSE,
          searching = FALSE,
          info = FALSE,
          scrollX = TRUE,
          dom = "t"
        )
      )
  })
  
  # PRODUCTION ENVELOPE PLOT #######################
  output$fba_production_envelope_plot <- renderPlotly({
    req(production_envelope)
    
    plot_ly(data = production_envelope) %>%
      add_trace(
        x = ~Growth,
        y = ~flux_maximum,
        name = "Max PHB Flux",
        type = "scatter",
        mode = "lines",
        line = list(color = "#21918c", width = 2),
        fill = "tozeroy",
        fillcolor = "rgba(33, 145, 140, 0.2)"
      ) %>%
      add_trace(
        x = ~Growth,
        y = ~flux_minimum,
        name = "Min PHB Flux",
        type = "scatter",
        mode = "lines",
        line = list(color = "#440154", width = 1.5, dash = "dash")
      ) %>%
      layout(
        xaxis = list(title = "Growth Rate (1/h)"),
        yaxis = list(title = "PHB Flux (mmol/gDW/h)"),
        hovermode = "x unified",
        legend = list(orientation = "h", x = 0.05, y = 1.15)
      )
  })
  
  # PFBA FLUX DISTRIBUTION TABLE ###################
  output$fba_pfba_table <- renderDT({
    req(pfba_flux_distribution)
    
    pfba_flux_distribution %>%
      transmute(
        `Reaction ID` = reaction_id,
        `Reaction Name` = reaction_name,
        `Flux (mmol/gDW/h)` = round(flux, 4)
      ) %>%
      datatable(
        rownames = FALSE,
        selection = "single",
        options = list(
          pageLength = 15,
          scrollX = TRUE
        )
      )
  })
  
}

###################################################################################################