###################################################################################################
# SERVER: Environmental Screening
###################################################################################################

page_environmental_screening_server <- function(input, output, session,
                                                  medium_screen_results,
                                                  minimal_medium_exchanges,
                                                  aeration_carbon_results) {

  # -------------------------------------------------------------------------------------------
  # TAB 1: MEDIUM SCREENING
  # -------------------------------------------------------------------------------------------

  output$env_medium_comparison_plot <- renderPlotly({
    req(medium_screen_results)

    df <- medium_screen_results %>%
      mutate(medium = factor(medium, levels = medium))

    plot_ly(data = df, x = ~medium) %>%
      add_trace(
        y = ~max_growth_rate_1_h,
        name = "Max Growth Rate (1/h)",
        type = "bar",
        marker = list(color = "#21918c"),
        yaxis = "y"
      ) %>%
      add_trace(
        y = ~growth_coupled_PHB_flux,
        name = "Growth-Coupled PHB Flux (mmol/gDW/h)",
        type = "bar",
        marker = list(color = "#440154"),
        yaxis = "y2"
      ) %>%
      layout(
        barmode = "group",
        xaxis  = list(title = "", tickangle = -20),
        yaxis  = list(title = "Growth Rate (1/h)"),
        yaxis2 = list(title = "PHB Flux (mmol/gDW/h)", overlaying = "y", side = "right"),
        legend = list(orientation = "h", x = 0.05, y = 1.15),
        margin = list(b = 100)
      )
  })

  output$env_medium_table <- renderDT({
    req(medium_screen_results)

    medium_screen_results %>%
      transmute(
        `Medium`                      = medium,
        `Glucose Uptake Bound`        = glc_uptake_bound,
        `NH4 Uptake Bound`            = nh4_uptake_bound,
        `O2 Uptake Bound`             = o2_uptake_bound,
        `Max Growth Rate (1/h)`       = round(max_growth_rate_1_h, 4),
        `Growth-Coupled Growth (1/h)` = round(growth_coupled_growth_1_h, 4),
        `Growth-Coupled PHB Flux`     = round(growth_coupled_PHB_flux, 4)
      ) %>%
      datatable(
        rownames  = FALSE,
        selection = "single",
        options   = list(paging = FALSE, searching = FALSE, info = FALSE, scrollX = TRUE, dom = "t")
      )
  })

  output$env_minimal_medium_table <- renderDT({
    req(minimal_medium_exchanges)

    minimal_medium_exchanges %>%
      transmute(
        `Exchange Reaction`        = X,
        `Uptake Flux (mmol/gDW/h)` = round(uptake_flux, 4)
      ) %>%
      arrange(desc(`Uptake Flux (mmol/gDW/h)`)) %>%
      datatable(
        rownames  = FALSE,
        selection = "single",
        options   = list(paging = FALSE, scrollX = TRUE, scrollY = "300px", dom = "t")
      )
  })

  # -------------------------------------------------------------------------------------------
  # TAB 2: AERATION x CARBON SOURCE
  # -------------------------------------------------------------------------------------------

  aeration_levels <- c("Anaerobic", "Microaerobic", "Aerobic", "Hyperaerobic")

  aeration_df <- reactive({
    req(aeration_carbon_results)

    aeration_carbon_results %>%
      mutate(aeration = factor(aeration, levels = aeration_levels))
  })

  # Reshapes the long (carbon_source x aeration) table into a matrix suitable for a plotly heatmap
  build_env_heatmap <- function(df, value_col, colorscale, title) {
    wide <- df %>%
      select(carbon_source, aeration, value = all_of(value_col)) %>%
      pivot_wider(names_from = carbon_source, values_from = value) %>%
      arrange(aeration)

    x_vals <- sort(setdiff(names(wide), "aeration"))
    y_vals <- as.character(wide$aeration)
    z_mat  <- as.matrix(wide[, x_vals])

    plot_ly(
      x = x_vals, y = y_vals, z = z_mat,
      type = "heatmap",
      colorscale = colorscale,
      colorbar = list(title = title)
    ) %>%
      layout(
        xaxis = list(title = "Carbon Source", tickangle = -45),
        yaxis = list(title = "Aeration Level")
      )
  }

  output$env_aeration_growth_heatmap <- renderPlotly({
    df <- aeration_df()
    shiny::validate(shiny::need(nrow(df) > 0, "No aeration / carbon-source data available."))

    build_env_heatmap(df, "max_growth_rate_1_h", "Viridis", "Growth Rate<br>(1/h)")
  })

  output$env_aeration_phb_heatmap <- renderPlotly({
    df <- aeration_df()
    shiny::validate(shiny::need(nrow(df) > 0, "No aeration / carbon-source data available."))

    build_env_heatmap(df, "growth_coupled_PHB_flux", "Plasma", "PHB Flux<br>(mmol/gDW/h)")
  })

  output$env_carbon_source_selector <- renderUI({
    req(aeration_carbon_results)

    choices <- sort(unique(aeration_carbon_results$carbon_source))

    selectInput(
      inputId  = "env_carbon_source_filter",
      label    = "Carbon Source:",
      choices  = choices,
      selected = if ("Glucose" %in% choices) "Glucose" else choices[1]
    )
  })

  output$env_aeration_table <- renderDT({
    req(input$env_carbon_source_filter)

    df <- aeration_df() %>%
      filter(carbon_source == input$env_carbon_source_filter) %>%
      arrange(aeration)

    shiny::validate(
      shiny::need(nrow(df) > 0, "No data available for this carbon source.")
    )

    df %>%
      transmute(
        `Carbon Source`           = carbon_source,
        `Aeration`                = as.character(aeration),
        `O2 Uptake Bound`         = o2_uptake_bound,
        `Max Growth Rate (1/h)`   = round(max_growth_rate_1_h, 4),
        `Growth-Coupled PHB Flux` = round(growth_coupled_PHB_flux, 4)
      ) %>%
      datatable(
        rownames  = FALSE,
        selection = "single",
        options   = list(paging = FALSE, searching = FALSE, info = FALSE, scrollX = TRUE, dom = "t")
      )
  })

}

###################################################################################################
