###################################################################################################
# SERVER: Genome Annotation
###################################################################################################

page_genome_server <- function(input, output, session, annotated_genome) {
  
  # SELECTED GENE ID ################################
  selected_gene_id <- reactive({
    req(input$genome_overview_table_rows_selected)
    idx <- input$genome_overview_table_rows_selected
    
    df <- annotated_genome$overview
    df$ID[idx]
  })
  
  # OVERVIEW TABLE ##################################
  output$genome_overview_table <- renderDT({
    req(annotated_genome$overview)
    
    datatable(
      annotated_genome$overview,
      rownames = FALSE,
      selection = "single",
      options = list(
        pageLength = 15,
        scrollX = TRUE
      )
    )
  })
  
  # DATABASE CROSS-REFERENCES TABLE #################
  output$genome_databases_table <- renderDT({
    selected_id <- selected_gene_id()
    
    shiny::validate(
      shiny::need(selected_id, "Please select a gene from the overview table.")
    )
    
    db_df <- annotated_genome$databases
    
    db_data <- db_df %>%
      filter(.data$ID == selected_id)
    
    shiny::validate(
      shiny::need(nrow(db_data) > 0, "No database cross-references available for this gene.")
    )
    
    datatable(
      db_data,
      rownames = FALSE,
      selection = "single",
      options = list(
        paging = FALSE,
        scrollX = TRUE,
        scrollY = "300px",
        dom = "t"
      )
    )
  })
  
  # BLAST RESULTS TABLE #############################
  output$genome_blast_table <- renderDT({
    selected_id <- selected_gene_id()
    
    shiny::validate(
      shiny::need(selected_id, "Please select a gene from the overview table.")
    )
    
    blast_df <- annotated_genome$blast
    
    blast_data <- blast_df %>%
      filter(.data$ID == selected_id)
    
    shiny::validate(
      shiny::need(nrow(blast_data) > 0, "No BLAST results available for this gene.")
    )
    
    datatable(
      blast_data,
      rownames = FALSE,
      selection = "single",
      options = list(
        paging = FALSE,
        scrollX = TRUE,
        scrollY = "300px",
        dom = "t"
      )
    )
  })
  
}

###################################################################################################