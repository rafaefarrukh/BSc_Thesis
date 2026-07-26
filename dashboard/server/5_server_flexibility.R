###################################################################################################
# SERVER: PHB Pathway Flexibility
###################################################################################################

page_flexibility_server <- function(input, output, session, pathway_flexibility) {
  
  # FILTERED DATASET ################################
  filtered_flexibility <- reactive({
    req(pathway_flexibility)
    
    df <- pathway_flexibility %>%
      filter(
        distance_from_PHB <= input$flexibility_max_distance,
        flux_range >= input$flexibility_min_range
      )
    
    if (nzchar(trimws(input$flexibility_search))) {
      pattern <- trimws(input$flexibility_search)
      df <- df %>%
        filter(
          str_detect(reaction_id, regex(pattern, ignore_case = TRUE)) |
            str_detect(reaction_name, regex(pattern, ignore_case = TRUE)) |
            str_detect(subsystem, regex(pattern, ignore_case = TRUE))
        )
    }
    
    df
  })
  
  # BAR CHART PLOT ##################################
  output$flexibility_bar_plot <- renderPlotly({
    df <- filtered_flexibility()
    
    shiny::validate(
      shiny::need(nrow(df) > 0, "No reactions match the selected filter criteria.")
    )
    
    plot_ly(
      data = df,
      x = ~reaction_id,
      y = ~flux_range,
      type = "bar",
      marker = list(
        color = ~distance_from_PHB,
        colorscale = "viridis",
        showscale = TRUE,
        colorbar = list(title = "Distance from<br>PHB Synthesis")
      ),
      hovertext = ~paste0(
        "<b>", reaction_id, "</b><br>",
        "Name: ", ifelse(is.na(reaction_name), "N/A", reaction_name), "<br>",
        "Subsystem: ", ifelse(is.na(subsystem), "N/A", subsystem), "<br>",
        "Distance: ", distance_from_PHB, "<br>",
        "Min Flux: ", round(minimum, 4), "<br>",
        "Max Flux: ", round(maximum, 4), "<br>",
        "Flux Range: ", round(flux_range, 4)
      ),
      hoverinfo = "text"
    ) %>%
      layout(
        xaxis = list(
          title = "Reaction ID",
          tickangle = -90,
          categoryorder = "trace"
        ),
        yaxis = list(
          title = "Feasible Flux Range (max - min, mmol/gDW/h)"
        ),
        margin = list(b = 100)
      )
  })
  
  # DATA TABLE ######################################
  output$flexibility_table <- renderDT({
    df <- filtered_flexibility()
    
    shiny::validate(
      shiny::need(nrow(df) > 0, "No reactions match the selected filter criteria.")
    )
    
    df %>%
      transmute(
        `Reaction ID`   = reaction_id,
        `Reaction Name` = ifelse(is.na(reaction_name), "—", reaction_name),
        `Subsystem`     = ifelse(is.na(subsystem), "—", subsystem),
        `Distance`      = distance_from_PHB,
        `Min Flux`      = round(minimum, 4),
        `Max Flux`      = round(maximum, 4),
        `Flux Range`    = round(flux_range, 4)
      ) %>%
      datatable(
        rownames = FALSE,
        selection = "single",
        options = list(
          pageLength = 15,
          scrollX    = TRUE
        )
      )
  })
  
}

###################################################################################################