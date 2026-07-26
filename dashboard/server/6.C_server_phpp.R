###################################################################################################
# SERVER: Phenotype Phase Plane (PhPP)
###################################################################################################

page_phpp_server <- function(input, output, session,
                              phpp_glucose_grid, phpp_glucose_transitions,
                              phpp_maltose_grid, phpp_maltose_transitions) {

  phpp_selected <- reactive({
    if (identical(input$phpp_carbon_source, "Maltose")) {
      list(
        grid    = phpp_maltose_grid,
        trans   = phpp_maltose_transitions,
        xcol    = "EX_malt_e",
        xlabel  = "Maltose Uptake (mmol/gDW/h)"
      )
    } else {
      list(
        grid    = phpp_glucose_grid,
        trans   = phpp_glucose_transitions,
        xcol    = "EX_glc__D_e",
        xlabel  = "Glucose Uptake (mmol/gDW/h)"
      )
    }
  })

  output$phpp_contour_plot <- renderPlotly({
    sel   <- phpp_selected()
    grid  <- sel$grid
    trans <- sel$trans

    shiny::validate(shiny::need(nrow(grid) > 0, "No PhPP grid data available for this carbon source."))

    # Reshape the long-format (carbon uptake, O2 uptake, growth rate) grid into a matrix.
    # NOTE: the grid's growth-rate column is named `flux_maximum` in the underlying CSV even
    # though it represents the maximum specific growth rate at each grid coordinate (per the
    # PhPP methodology), not a PHB flux.
    wide <- grid %>%
      select(xv = all_of(sel$xcol), yv = EX_o2_e, value = flux_maximum) %>%
      pivot_wider(names_from = xv, values_from = value) %>%
      arrange(yv)

    value_cols <- setdiff(names(wide), "yv")
    value_cols <- value_cols[order(as.numeric(value_cols))]

    z_mat  <- as.matrix(wide[, value_cols])
    x_vals <- as.numeric(value_cols)
    y_vals <- wide$yv

    plot_ly() %>%
      add_trace(
        x = x_vals, y = y_vals, z = z_mat,
        type = "contour",
        colorscale = "Viridis",
        contours = list(showlabels = TRUE, coloring = "heatmap"),
        colorbar = list(title = "Growth Rate<br>(1/h)"),
        name = "Growth rate surface"
      ) %>%
      add_trace(
        x = trans[[sel$xcol]],
        y = trans$EX_o2_e,
        type = "scatter",
        mode = "markers",
        marker = list(color = "red", size = 7, symbol = "x"),
        name = "Candidate phase transition"
      ) %>%
      layout(
        xaxis  = list(title = sel$xlabel),
        yaxis  = list(title = "Oxygen Uptake (mmol/gDW/h)"),
        legend = list(orientation = "h", x = 0.05, y = 1.1)
      )
  })

  output$phpp_region_table <- renderDT({
    sel   <- phpp_selected()
    trans <- sel$trans

    shiny::validate(shiny::need(nrow(trans) > 0, "No candidate transition points available."))

    trans %>%
      transmute(
        `Carbon Uptake`      = round(.data[[sel$xcol]], 3),
        `O2 Uptake`          = round(EX_o2_e, 3),
        `Growth Rate (1/h)`  = round(growth_rate, 4),
        `Curvature`          = signif(curvature, 4)
      ) %>%
      arrange(desc(Curvature)) %>%
      datatable(
        rownames  = FALSE,
        selection = "single",
        options   = list(pageLength = 15, scrollX = TRUE)
      )
  })

}

###################################################################################################
