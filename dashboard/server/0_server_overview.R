###################################################################################################
# SERVER: Overview
###################################################################################################

page_overview_server <- function(input, output, session, annotated_genome, gem_network) {

  # GENOME / NETWORK SUMMARY STATS ##################
  # Static snapshot stats derived from the same data objects used on the Genome and GEM
  # pages, so this page never needs its own separate data source for these three boxes.

  output$overview_genes_out <- renderText({
    format(nrow(annotated_genome$overview), big.mark = ",")
  })

  output$overview_reactions_out <- renderText({
    format(nrow(gem_network$reaction_nodes), big.mark = ",")
  })

  output$overview_metabolites_out <- renderText({
    format(nrow(gem_network$metabolite_nodes), big.mark = ",")
  })

  # COST / YIELD SUMMARY ############################
  # Mirrors the default-file loading pattern used in page_cost_server(), but computed once
  # here as a read-only snapshot across every medium x recovery-method combination, since
  # the Overview page doesn't need to react to edits made on the Cost Analysis page.

  overview_cost_summary <- reactive({
    medium_df   <- tryCatch(read_csv("data/8_media.csv",      show_col_types = FALSE), error = function(e) NULL)
    recovery_df <- tryCatch(read_csv("data/8_downstream.csv", show_col_types = FALSE), error = function(e) NULL)
    yield_df    <- tryCatch(read_csv("data/8_phb.csv",        show_col_types = FALSE), error = function(e) NULL)

    if (is.null(medium_df) || is.null(recovery_df) || is.null(yield_df)) {
      return(NULL)
    }

    yield_col <- if ("PHB" %in% colnames(yield_df)) "PHB" else "PHB (g/1L volume)"
    cost_col  <- if ("Cost" %in% colnames(recovery_df)) "Cost" else "Cost (PKR/kg PHB)"

    media_names <- setdiff(colnames(medium_df), c("Component", "Cost"))

    per_medium <- lapply(media_names, function(m) {
      amounts <- suppressWarnings(as.numeric(medium_df[[m]]))
      costs   <- suppressWarnings(as.numeric(medium_df[["Cost"]]))
      medium_cost_L <- sum(amounts * costs / 1000, na.rm = TRUE)

      y_row <- yield_df[yield_df$`Medium Name` == m, ]
      y_val <- if (nrow(y_row) > 0) suppressWarnings(as.numeric(y_row[[yield_col]][1])) else NA_real_
      if (is.na(y_val) || y_val <= 0) y_val <- NA_real_

      medium_cost_kg <- if (!is.na(y_val)) (medium_cost_L / y_val) * 1000 else NA_real_

      recovery_costs_kg <- suppressWarnings(as.numeric(recovery_df[[cost_col]])) *
        suppressWarnings(as.numeric(recovery_df[["Scaling Factor"]]))
      recovery_cost_kg <- if (!is.na(y_val)) (min(recovery_costs_kg, na.rm = TRUE) / y_val) * 1000 else NA_real_

      tibble(
        medium    = m,
        yield_gL  = y_val,
        total_cost = medium_cost_kg + recovery_cost_kg
      )
    })

    bind_rows(per_medium) %>% filter(!is.na(total_cost))
  })

  output$overview_best_yield_out <- renderText({
    df <- overview_cost_summary()
    if (is.null(df) || nrow(df) == 0) return("N/A")

    best <- df %>% arrange(desc(yield_gL)) %>% slice(1)
    paste0(round(best$yield_gL, 2), " g/L (", best$medium, ")")
  })

  output$overview_best_cost_out <- renderText({
    df <- overview_cost_summary()
    if (is.null(df) || nrow(df) == 0) return("N/A")

    best <- df %>% arrange(total_cost) %>% slice(1)
    paste0(format(round(best$total_cost, 0), big.mark = ","), " PKR/kg (", best$medium, ")")
  })

}

###################################################################################################
