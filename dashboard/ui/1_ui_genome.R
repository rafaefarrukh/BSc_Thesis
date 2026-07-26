###################################################################################################
# UI: Genome Annotation
###################################################################################################

page_genome <- function() {
  layout_columns(
    col_widths = c(7, 5),
    
    card(
      full_screen = TRUE,
      card_header("Genome Annotation Overview"),
      DTOutput("genome_overview_table")
    ),
    
    layout_column_wrap(
      width = 1,
      card(
        full_screen = TRUE,
        card_header("Database Cross-References"),
        DTOutput("genome_databases_table")
      ),
      card(
        full_screen = TRUE,
        card_header("BLAST Results"),
        DTOutput("genome_blast_table")
      )
    )
  )
}

###################################################################################################
