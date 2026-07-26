###################################################################################################
# UI: Overview
###################################################################################################

page_overview <- function() {
  tagList(

    # HEADER / DESCRIPTION ##########################
    card(
      height = 1/4,
      card_header("In-Silico Optimization of Polyhydroxyalkanoate (PHB) Production"),
      markdown("
      **Rossellomorea marisflavi ARS23** — a halophilic, PHB-producing soil isolate — was
      characterized experimentally and used to reconstruct a curated genome-scale metabolic
      model (GEM). The model was analyzed to understand PHB biosynthesis and to identify
      medium, environmental, and strain-design strategies for improving yield and reducing cost.

      Use the tabs above to explore each stage of the pipeline: genome annotation, the
      predicted phaC enzyme structure, the metabolic network, flux balance analyses,
      sensitivity/perturbation studies, system diagnostics, strain design, and the final
      cost analysis.
      ")
    ),

    # KEY METRICS ####################################
    layout_column_wrap(
      height = 1/4,
      value_box(
        title = "Genes Annotated",
        value = textOutput("overview_genes_out"),
        showcase = icon("dna"),
        theme = "primary"
      ),
      value_box(
        title = "Reactions in Network",
        value = textOutput("overview_reactions_out"),
        showcase = icon("diagram-project"),
        theme = "secondary"
      ),
      value_box(
        title = "Metabolites in Network",
        value = textOutput("overview_metabolites_out"),
        showcase = icon("circle-nodes"),
        theme = "secondary"
      ),
      value_box(
        title = "Best PHB Titer",
        value = textOutput("overview_best_yield_out"),
        showcase = icon("flask-vial"),
        theme = "success"
      ),
      value_box(
        title = "Lowest Production Cost",
        value = textOutput("overview_best_cost_out"),
        showcase = icon("coins"),
        theme = "success"
      )
    ),

    # PIPELINE SUMMARY ################################
    card(
      height = 2/4,
      card_header("Methodological Pipeline"),
      layout_columns(
        col_widths = c(3, 3, 3, 3),

        card(
          class = "bg-light",
          card_body(
            tags$h6("1. Screening & Characterization"),
            tags$p(class = "small text-muted",
              "Sudan Black B / Nile Blue A screening, FTIR confirmation, and RSM-based ",
              "optimization of growth conditions (temperature, pH, NaCl)."
            )
          )
        ),
        card(
          class = "bg-light",
          card_body(
            tags$h6("2. GEM Construction"),
            tags$p(class = "small text-muted",
              "Genome annotation (Bakta), phaC structural validation (ColabFold), and a ",
              "bottom-up (KBase) + top-down (CarveMe) draft reconstruction."
            )
          )
        ),
        card(
          class = "bg-light",
          card_body(
            tags$h6("3. Curation & Validation"),
            tags$p(class = "small text-muted",
              "Manual curation of PHB pathway, annotations, and unbounded flux, validated ",
              "at each stage with MEMOTE."
            )
          )
        ),
        card(
          class = "bg-light",
          card_body(
            tags$h6("4. Analyses & Cost"),
            tags$p(class = "small text-muted",
              "FBA/pFBA, sensitivity and environmental perturbations, system diagnostics, ",
              "in silico strain design, and downstream cost analysis."
            )
          )
        )
      )
    )
  )
}

###################################################################################################
