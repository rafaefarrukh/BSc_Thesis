###################################################################################################
# SERVER: Dynamic FBA (dFBA)
###################################################################################################

page_dfba_server <- function(input, output, session,
                              dfba_glucose_standard, dfba_glucose_nlimited, dfba_glucose_twostage,
                              dfba_maltose_standard, dfba_maltose_nlimited, dfba_maltose_twostage) {

  dfba_data <- reactive({
    req(input$dfba_carbon, input$dfba_condition)

    carbon <- input$dfba_carbon
    cond   <- input$dfba_condition

    df <- switch(
      paste(carbon, cond, sep = "_"),
      "glucose_standard" = dfba_glucose_standard,
      "glucose_nlimited" = dfba_glucose_nlimited,
      "glucose_twostage" = dfba_glucose_twostage,
      "maltose_standard" = dfba_maltose_standard,
      "maltose_nlimited" = dfba_maltose_nlimited,
      "maltose_twostage" = dfba_maltose_twostage
    )

    substrate_col  <- if (carbon == "glucose") "glucose_mM" else "maltose_mM"
    substrate_name <- if (carbon == "glucose") "Glucose (mM)" else "Maltose (mM)"

    df %>%
      rename(substrate_mM = all_of(substrate_col)) %>%
      mutate(substrate_name = substrate_name)
  })

  # TIME-SERIES PLOT (2x2 grid: biomass, substrate, ammonium, PHB) ############################

  output$dfba_timeseries_plot <- renderPlotly({
    df <- dfba_data()

    shiny::validate(shiny::need(nrow(df) > 0, "No dFBA data available for this selection."))

    # Standard / nitrogen-limited runs diverge to non-physiological magnitudes past ~80h
    # (explicit-Euler numerical artifact, see warning card); truncate the plot accordingly.
    is_twostage <- identical(input$dfba_condition, "twostage")
    plot_until  <- if (is_twostage) max(df$time_h) else min(80, max(df$time_h))
    df_plot     <- df %>% filter(time_h <= plot_until)

    substrate_label <- unique(df_plot$substrate_name)[1]

    p_biomass <- plot_ly(df_plot, x = ~time_h, y = ~biomass_gL, type = "scatter", mode = "lines",
                          line = list(color = "#21918c"), name = "Biomass (g/L)")
    p_substrate <- plot_ly(df_plot, x = ~time_h, y = ~substrate_mM, type = "scatter", mode = "lines",
                            line = list(color = "#3b528b"), name = substrate_label)
    p_ammonium <- plot_ly(df_plot, x = ~time_h, y = ~ammonium_mM, type = "scatter", mode = "lines",
                           line = list(color = "#5ec962"), name = "Ammonium (mM)")
    p_phb <- plot_ly(df_plot, x = ~time_h, y = ~PHB_mM, type = "scatter", mode = "lines",
                      line = list(color = "#440154"), name = "PHB (mM)")

    subplot(p_biomass, p_substrate, p_ammonium, p_phb, nrows = 2, shareX = TRUE, titleY = TRUE, margin = 0.07) %>%
      layout(
        showlegend = FALSE,
        annotations = list(
          list(x = 0.18, y = 1.00, xref = "paper", yref = "paper", xanchor = "center", showarrow = FALSE, text = "<b>Biomass</b>"),
          list(x = 0.82, y = 1.00, xref = "paper", yref = "paper", xanchor = "center", showarrow = FALSE, text = paste0("<b>", substrate_label, "</b>")),
          list(x = 0.18, y = 0.44, xref = "paper", yref = "paper", xanchor = "center", showarrow = FALSE, text = "<b>Ammonium</b>"),
          list(x = 0.82, y = 0.44, xref = "paper", yref = "paper", xanchor = "center", showarrow = FALSE, text = "<b>PHB</b>")
        )
      )
  })

  # DATA TABLE #################################################################################

  output$dfba_table <- renderDT({
    df <- dfba_data()

    shiny::validate(shiny::need(nrow(df) > 0, "No dFBA data available for this selection."))

    substrate_label <- unique(df$substrate_name)[1]

    df %>%
      transmute(
        `Time (h)`      = round(time_h, 2),
        `Biomass (g/L)` = signif(biomass_gL, 5),
        !!substrate_label := signif(substrate_mM, 5),
        `Ammonium (mM)` = signif(ammonium_mM, 5),
        `PHB (mM)`      = signif(PHB_mM, 5)
      ) %>%
      datatable(
        rownames  = FALSE,
        selection = "single",
        options   = list(pageLength = 15, scrollX = TRUE)
      )
  })

}

###################################################################################################
