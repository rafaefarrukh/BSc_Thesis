###################################################################################################
# SERVER: Multi-Fraction & Loopless FVA
###################################################################################################

page_multifva_server <- function(input, output, session,
                                  multifva_100_standard, multifva_100_loopless,
                                  multifva_90_standard,  multifva_90_loopless,
                                  multifva_50_standard,  multifva_50_loopless,
                                  pathway_flexibility) {

  # Reaction metadata (name / subsystem / distance) is not present in the multifva_* tables
  # themselves (they only carry reaction id + min/max flux) - it's sourced from
  # pathway_flexibility, which covers the same 150-reaction PHB pathway neighborhood.
  reaction_meta <- pathway_flexibility %>%
    select(reaction_id, reaction_name, subsystem, distance_from_PHB) %>%
    distinct(reaction_id, .keep_all = TRUE)

  get_fraction_df <- function(fraction, method) {
    df <- switch(
      paste(fraction, method, sep = "_"),
      "100_standard" = multifva_100_standard,
      "100_loopless" = multifva_100_loopless,
      "90_standard"  = multifva_90_standard,
      "90_loopless"  = multifva_90_loopless,
      "50_standard"  = multifva_50_standard,
      "50_loopless"  = multifva_50_loopless
    )

    df %>%
      rename(reaction_id = X) %>%
      mutate(
        fraction   = paste0(fraction, "% of PHB optimum"),
        flux_range = maximum - minimum
      )
  }

  multifva_combined <- reactive({
    req(input$multifva_fractions, input$multifva_method)

    combined <- bind_rows(
      lapply(input$multifva_fractions, get_fraction_df, method = input$multifva_method)
    )

    combined <- combined %>%
      left_join(reaction_meta, by = "reaction_id") %>%
      filter(
        distance_from_PHB <= input$multifva_max_distance,
        flux_range >= input$multifva_min_range
      )

    if (nzchar(trimws(input$multifva_search))) {
      pattern <- trimws(input$multifva_search)
      combined <- combined %>%
        filter(
          str_detect(reaction_id, regex(pattern, ignore_case = TRUE)) |
            str_detect(coalesce(reaction_name, ""), regex(pattern, ignore_case = TRUE)) |
            str_detect(coalesce(as.character(subsystem), ""), regex(pattern, ignore_case = TRUE))
        )
    }

    combined
  })

  output$multifva_bar_plot <- renderPlotly({
    df <- multifva_combined()

    shiny::validate(shiny::need(nrow(df) > 0, "No reactions match the selected filter criteria."))

    plot_ly(
      data = df,
      x = ~reaction_id,
      y = ~flux_range,
      color = ~fraction,
      type = "bar",
      hovertext = ~paste0(
        "<b>", reaction_id, "</b><br>",
        "Fraction: ", fraction, "<br>",
        "Min: ", round(minimum, 4), "<br>",
        "Max: ", round(maximum, 4), "<br>",
        "Range: ", round(flux_range, 4)
      ),
      hoverinfo = "text"
    ) %>%
      layout(
        barmode = "group",
        xaxis   = list(title = "Reaction ID", tickangle = -90, categoryorder = "trace"),
        yaxis   = list(title = "Flux Range (max - min, mmol/gDW/h)"),
        margin  = list(b = 120),
        legend  = list(orientation = "h", x = 0.05, y = 1.15)
      )
  })

  output$multifva_table <- renderDT({
    df <- multifva_combined()

    shiny::validate(shiny::need(nrow(df) > 0, "No reactions match the selected filter criteria."))

    df %>%
      transmute(
        `Reaction ID`   = reaction_id,
        `Reaction Name` = ifelse(is.na(reaction_name) | reaction_name == "", "\u2014", reaction_name),
        `Subsystem`     = ifelse(is.na(subsystem), "\u2014", as.character(subsystem)),
        `Distance`      = distance_from_PHB,
        `Fraction`      = fraction,
        `Min Flux`      = round(minimum, 4),
        `Max Flux`      = round(maximum, 4),
        `Flux Range`    = round(flux_range, 4)
      ) %>%
      arrange(`Reaction ID`, `Fraction`) %>%
      datatable(
        rownames  = FALSE,
        selection = "single",
        options   = list(pageLength = 15, scrollX = TRUE)
      )
  })

}

###################################################################################################
