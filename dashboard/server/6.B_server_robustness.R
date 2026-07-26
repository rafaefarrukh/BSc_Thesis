###################################################################################################
# SERVER: Robustness (1D & 2D Sensitivity Analyses)
###################################################################################################

page_robustness_server <- function(input, output, session,
                                    robust_1d_glucose, robust_1d_maltose,
                                    robust_1d_nitrogen, robust_1d_oxygen,
                                    robust_2d_glc_nh4, robust_2d_malt_nh4, robust_2d_glc_o2) {

  # -------------------------------------------------------------------------------------------
  # 1D NUTRIENT SWEEP
  # -------------------------------------------------------------------------------------------

  robust_1d_selected <- reactive({
    switch(
      input$robust_nutrient,
      "O2"   = list(df = robust_1d_oxygen,   label = "O2 Uptake Bound (mmol/gDW/h)"),
      "Glc"  = list(df = robust_1d_glucose,  label = "Glucose Uptake Bound (mmol/gDW/h)"),
      "Malt" = list(df = robust_1d_maltose,  label = "Maltose Uptake Bound (mmol/gDW/h)"),
      "NH4"  = list(df = robust_1d_nitrogen, label = "Ammonia (NH4) Uptake Bound (mmol/gDW/h)")
    )
  })

  output$robust_line_plot <- renderPlotly({
    sel <- robust_1d_selected()
    df  <- sel$df

    shiny::validate(shiny::need(nrow(df) > 0, "No sweep data available for this nutrient."))

    plot_ly(data = df, x = ~uptake_bound) %>%
      add_trace(
        y = ~growth_rate_1_h,
        name = "Max Growth Rate (1/h)",
        type = "scatter",
        mode = "lines+markers",
        line = list(color = "#21918c", width = 2),
        marker = list(color = "#21918c"),
        yaxis = "y"
      ) %>%
      add_trace(
        y = ~PHB_flux,
        name = "Growth-Coupled PHB Flux (mmol/gDW/h)",
        type = "scatter",
        mode = "lines+markers",
        line = list(color = "#440154", width = 2, dash = "dash"),
        marker = list(color = "#440154"),
        yaxis = "y2"
      ) %>%
      layout(
        xaxis  = list(title = sel$label),
        yaxis  = list(title = "Growth Rate (1/h)"),
        yaxis2 = list(title = "PHB Flux (mmol/gDW/h)", overlaying = "y", side = "right"),
        hovermode = "x unified",
        legend = list(orientation = "h", x = 0.05, y = 1.15)
      )
  })

  # -------------------------------------------------------------------------------------------
  # 2D NUTRIENT GRID HEATMAP
  # -------------------------------------------------------------------------------------------

  # Reshapes a long-format (x, y, value) grid into a matrix for a plotly heatmap. Column names
  # are sorted numerically (not by round-tripping through as.character()) to avoid mismatches
  # between the pivoted column-name strings and re-stringified numeric values.
  build_robust_heatmap <- function(df, xcol, ycol, value_col, colorscale, title, xlab, ylab) {
    wide <- df %>%
      select(xv = all_of(xcol), yv = all_of(ycol), value = all_of(value_col)) %>%
      pivot_wider(names_from = xv, values_from = value) %>%
      arrange(yv)

    value_cols <- setdiff(names(wide), "yv")
    value_cols <- value_cols[order(as.numeric(value_cols))]

    z_mat  <- as.matrix(wide[, value_cols])
    x_vals <- as.numeric(value_cols)
    y_vals <- wide$yv

    plot_ly(
      x = x_vals, y = y_vals, z = z_mat,
      type = "heatmap",
      colorscale = colorscale,
      colorbar = list(title = title)
    ) %>%
      layout(xaxis = list(title = xlab), yaxis = list(title = ylab))
  }

  output$robust_heatmap_plot <- renderPlotly({
    sel <- switch(
      input$robust_pair,
      "glc_o2"   = list(df = robust_2d_glc_o2,   xcol = "glucose", ycol = "oxygen",
                         xlab = "Glucose Uptake Bound", ylab = "Oxygen Uptake Bound"),
      "glc_nh4"  = list(df = robust_2d_glc_nh4,  xcol = "glucose", ycol = "nitrogen",
                         xlab = "Glucose Uptake Bound", ylab = "Ammonia (NH4) Uptake Bound"),
      "malt_nh4" = list(df = robust_2d_malt_nh4, xcol = "maltose", ycol = "nitrogen",
                         xlab = "Maltose Uptake Bound", ylab = "Ammonia (NH4) Uptake Bound")
    )

    shiny::validate(shiny::need(nrow(sel$df) > 0, "No grid data available for this nutrient pair."))

    growth_hm <- build_robust_heatmap(sel$df, sel$xcol, sel$ycol, "growth_rate_1_h",
                                       "Viridis", "Growth<br>(1/h)", sel$xlab, sel$ylab)
    phb_hm    <- build_robust_heatmap(sel$df, sel$xcol, sel$ycol, "PHB_flux",
                                       "Plasma", "PHB Flux<br>(mmol/gDW/h)", sel$xlab, "")

    subplot(growth_hm, phb_hm, nrows = 1, shareY = TRUE, titleX = TRUE, titleY = TRUE, margin = 0.08) %>%
      layout(
        annotations = list(
          list(x = 0.20, y = 1.08, xref = "paper", yref = "paper", showarrow = FALSE, text = "<b>Growth Rate</b>"),
          list(x = 0.80, y = 1.08, xref = "paper", yref = "paper", showarrow = FALSE, text = "<b>Growth-Coupled PHB Flux</b>")
        )
      )
  })

}

###################################################################################################
