###################################################################################################
# SERVER: Pathway Bottlenecks
###################################################################################################

page_bottlenecks_server <- function(input, output, session, bottleneck_results) {

  ranked <- reactive({
    req(bottleneck_results)

    bottleneck_results %>%
      arrange(desc(bottleneck_score))
  })

  output$diag_bottleneck_plot <- renderPlotly({
    df <- ranked() %>% slice_head(n = 20)

    shiny::validate(shiny::need(nrow(df) > 0, "No bottleneck data available."))

    plot_ly(
      data = df %>% arrange(bottleneck_score),
      x = ~bottleneck_score,
      y = ~factor(reaction_id, levels = reaction_id),
      type = "bar",
      orientation = "h",
      marker = list(
        color = ~distance_from_PHB,
        colorscale = "Viridis",
        showscale = TRUE,
        colorbar = list(title = "Distance from<br>PHB Synthesis")
      ),
      hovertext = ~paste0(
        "<b>", reaction_id, "</b><br>",
        reaction_name, "<br>",
        "Distance: ", distance_from_PHB, "<br>",
        "|Flux|: ", round(abs_flux, 3), "<br>",
        "Flux/Capacity: ", round(flux_to_capacity_ratio, 4), "<br>",
        "Bottleneck Score: ", round(bottleneck_score, 5)
      ),
      hoverinfo = "text"
    ) %>%
      layout(
        xaxis  = list(title = "Bottleneck Score"),
        yaxis  = list(title = ""),
        margin = list(l = 130)
      )
  })

  output$diag_bottleneck_table <- renderDT({
    df <- ranked()

    shiny::validate(shiny::need(nrow(df) > 0, "No bottleneck data available."))

    df %>%
      transmute(
        `Reaction ID`       = reaction_id,
        `Reaction Name`     = reaction_name,
        `Distance from PHB` = distance_from_PHB,
        `Flux`              = round(flux, 4),
        `|Flux|`            = round(abs_flux, 4),
        `Capacity Bound`    = capacity_bound,
        `Flux / Capacity`   = round(flux_to_capacity_ratio, 5),
        `Bottleneck Score`  = round(bottleneck_score, 5)
      ) %>%
      datatable(
        rownames  = FALSE,
        selection = "single",
        options   = list(pageLength = 15, scrollX = TRUE)
      )
  })

}

###################################################################################################
