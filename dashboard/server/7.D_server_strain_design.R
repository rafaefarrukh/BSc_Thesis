###################################################################################################
# SERVER: Strain Design
###################################################################################################

page_strain_design_server <- function(input, output, session,
                                       knockout_essentiality, knockout_single, knockout_double) {

  # -------------------------------------------------------------------------------------------
  # TAB 1: ESSENTIAL REACTIONS
  # -------------------------------------------------------------------------------------------

  output$strain_essential_summary_table <- renderDT({
    req(knockout_essentiality)

    knockout_essentiality %>%
      count(classification, name = "n") %>%
      mutate(
        classification = factor(
          classification,
          levels = c("essential", "critical", "non-essential"),
          labels = c("Essential", "Critical", "Non-essential")
        )
      ) %>%
      arrange(classification) %>%
      transmute(
        `Classification`  = classification,
        `Count`           = n,
        `% of Reactions`  = round(100 * n / sum(n), 1)
      ) %>%
      datatable(
        rownames  = FALSE,
        selection = "single",
        options   = list(paging = FALSE, searching = FALSE, info = FALSE, dom = "t")
      )
  })

  output$strain_essential_plot <- renderPlotly({
    req(knockout_essentiality)

    df <- knockout_essentiality %>%
      filter(classification %in% c("essential", "critical")) %>%
      arrange(growth_fraction_of_wt) %>%
      slice_head(n = 30)

    shiny::validate(shiny::need(nrow(df) > 0, "No essential/critical knockouts found."))

    plot_ly(
      data = df %>% arrange(desc(growth_fraction_of_wt)),
      x = ~growth_fraction_of_wt,
      y = ~factor(reaction_id, levels = reaction_id),
      color = ~classification,
      colors = c("essential" = "#c0392b", "critical" = "#e67e22"),
      type = "bar",
      orientation = "h",
      hovertext = ~paste0(
        "<b>", reaction_id, "</b><br>",
        ifelse(reaction_name == "", "", paste0(reaction_name, "<br>")),
        "Growth Fraction: ", round(growth_fraction_of_wt, 4), "<br>",
        "Classification: ", classification
      ),
      hoverinfo = "text"
    ) %>%
      layout(
        xaxis  = list(title = "Growth Fraction of Wild Type"),
        yaxis  = list(title = ""),
        margin = list(l = 150),
        legend = list(orientation = "h", x = 0.05, y = 1.08)
      )
  })

  # -------------------------------------------------------------------------------------------
  # TAB 2: GROWTH-COUPLED KNOCKOUTS
  # -------------------------------------------------------------------------------------------

  strain_knockout_df <- reactive({
    req(input$strain_knockout_type)

    df <- switch(
      input$strain_knockout_type,
      "single" = knockout_single,
      "double" = knockout_double
    )

    df %>% arrange(desc(min_PHB_at_max_growth))
  })

  output$strain_knockout_table <- renderDT({
    df <- strain_knockout_df()

    shiny::validate(shiny::need(nrow(df) > 0, "No knockout data available."))

    df %>%
      transmute(
        `Knockout(s)`                       = knockouts,
        `# Reactions Knocked Out`           = n_knockouts,
        `Max Growth Rate (1/h)`             = round(max_growth, 4),
        `Forced Min PHB Flux at Max Growth` = round(min_PHB_at_max_growth, 4)
      ) %>%
      datatable(
        rownames  = FALSE,
        selection = "single",
        options   = list(pageLength = 15, scrollX = TRUE)
      )
  })

}

###################################################################################################
