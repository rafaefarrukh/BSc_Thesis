###################################################################################################
# SERVER: GEM
###################################################################################################

page_gem_server <- function(input, output, session, gem_network) {

  # ACTIVE PRESET ##################################
  # "phb" and "complete" are implemented.
  #
  # gem_network$reaction_nodes already holds every reaction in the curated PHB
  # pathway distance table (i.e. the whole reachable network, not just a small
  # neighborhood), so "Complete Network" doesn't need its own dataset - it just
  # drives the max-distance slider up to the true maximum distance present in the
  # data, which removes the distance filter entirely.

  max_available_distance <- max(gem_network$reaction_nodes$distance, na.rm = TRUE)

  updateSliderInput(session, "gem_max_distance", max = max_available_distance)

  active_preset <- reactiveVal("phb")

  observeEvent(input$gem_preset_phb, {
    active_preset("phb")
    updateSliderInput(session, "gem_max_distance", value = 2)
  })

  observeEvent(input$gem_preset_complete, {
    active_preset("complete")
    updateSliderInput(session, "gem_max_distance", value = max_available_distance)
  })

  # FILTERED REACTION SET ##########################
  # Always scoped to reactions present in the PHB pathway distance table (never anything
  # outside of it), further filtered by the max-distance slider.

  filtered_reactions <- reactive({
    gem_network$reaction_nodes %>%
      filter(distance <= input$gem_max_distance)
  })

  filtered_edges <- reactive({
    rxn_ids <- filtered_reactions()$id

    edges <- gem_network$edges %>%
      filter(reaction_id %in% rxn_ids)

    if (isTRUE(input$gem_hide_currency)) {
      edges <- edges %>% filter(!metabolite_id %in% (gem_network$metabolite_nodes %>% filter(is_currency) %>% pull(id)))
    }

    edges
  })

  filtered_metabolites <- reactive({
    if (!isTRUE(input$gem_show_metabolites)) {
      return(gem_network$metabolite_nodes %>% filter(FALSE))
    }

    met_ids <- unique(filtered_edges()$metabolite_id)
    gem_network$metabolite_nodes %>% filter(id %in% met_ids)
  })

  # VISNETWORK NODE/EDGE TABLES ####################

  vis_nodes <- reactive({
    req(nrow(filtered_reactions()) > 0)

    dist_range <- gem_network$distance_range
    pal <- colorRampPalette(c("#fde725", "#21918c", "#440154"))(max(dist_range[2] - dist_range[1] + 1, 2))

    rxn_nodes <- filtered_reactions() %>%
      mutate(
        color_idx = distance - dist_range[1] + 1,
        color     = pal[pmin(pmax(color_idx, 1), length(pal))]
      ) %>%
      transmute(
        id       = node_id,
        label    = label,
        title    = paste0(
          "<b>", name, "</b><br>",
          "Distance from PHB: ", distance, "<br>",
          "Bounds: [", lower_bound, ", ", upper_bound, "]"
        ),
        group    = "reaction",
        color    = color,
        shape    = "box",
        font.size = 12,
        distance = distance
      )

    met_nodes <- filtered_metabolites() %>%
      transmute(
        id    = node_id,
        label = label,
        title = paste0("<b>", name, "</b><br>", "Compartment: ", compartment),
        group = "metabolite",
        color = "#d9d9d9",
        shape = "dot",
        font.size = 10,
        distance = NA_integer_
      )

    bind_rows(rxn_nodes, met_nodes)
  })

  vis_edges <- reactive({
    req(nrow(filtered_edges()) > 0)

    filtered_edges() %>%
      transmute(
        from   = from,
        to     = to,
        arrows = "to",
        color.color   = "#bbbbbb",
        color.opacity = 0.5
      )
  })

  # PHYSICS TOGGLE ###################################
  # Uses a proxy so toggling physics doesn't force a full re-render/re-stabilization
  # of the graph - it just freezes or resumes the existing layout in place.

  physics_enabled <- reactiveVal(TRUE)

  observeEvent(input$gem_toggle_physics, {
    physics_enabled(!physics_enabled())

    visNetworkProxy("gem_network_plot") %>%
      visPhysics(enabled = physics_enabled())

    updateActionButton(
      session,
      "gem_toggle_physics",
      label = if (physics_enabled()) "Freeze Layout" else "Resume Physics",
      icon  = icon(if (physics_enabled()) "lock" else "lock-open")
    )
  })

  # RENDER NETWORK ##################################

  output$gem_network_plot <- renderVisNetwork({
    nodes <- vis_nodes()
    edges <- vis_edges()

    visNetwork(nodes, edges) %>%
      visEdges(smooth = FALSE) %>%
      visOptions(
        highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
        nodesIdSelection = TRUE
      ) %>%
      visPhysics(
        solver = "forceAtlas2Based",
        enabled = isolate(physics_enabled()),
        stabilization = list(enabled = TRUE, iterations = 200),
        forceAtlas2Based = list(gravitationalConstant = -60, springLength = 90)
      ) %>%
      visInteraction(hover = TRUE, tooltipDelay = 100) %>%
      visEvents(selectNode = "function(nodes) {
        Shiny.setInputValue('gem_selected_node', nodes.nodes[0]);
      }")
  })

  # SELECTED NODE DETAILS ###########################

  output$gem_node_details <- renderUI({
    sel <- input$gem_selected_node
    if (is.null(sel) || length(sel) == 0) {
      return(div(class = "text-muted", "Click a node in the network to see details here."))
    }
    sel <- sel[1]

    if (startsWith(sel, "rxn_")) {
      rid <- sub("^rxn_", "", sel)
      row <- gem_network$reaction_nodes %>% filter(id == rid)

      if (nrow(row) == 0) {
        return(div(class = "text-muted", "Reaction not found."))
      }

      tagList(
        h5(row$name[1]),
        tags$table(
          class = "table table-sm",
          tags$tr(tags$td(strong("ID")),          tags$td(row$id[1])),
          tags$tr(tags$td(strong("Distance")),     tags$td(row$distance[1])),
          tags$tr(tags$td(strong("Bounds")),       tags$td(paste0("[", row$lower_bound[1], ", ", row$upper_bound[1], "]"))),
          tags$tr(tags$td(strong("GPR")),          tags$td(if (row$gene_reaction_rule[1] == "") em("none") else row$gene_reaction_rule[1])),
          tags$tr(tags$td(strong("Equation")),     tags$td(row$equation[1]))
        )
      )
    } else if (startsWith(sel, "met_")) {
      mid <- sub("^met_", "", sel)
      row <- gem_network$metabolite_nodes %>% filter(id == mid)

      if (nrow(row) == 0) {
        return(div(class = "text-muted", "Metabolite not found."))
      }

      tagList(
        h5(row$name[1]),
        tags$table(
          class = "table table-sm",
          tags$tr(tags$td(strong("ID")),          tags$td(row$id[1])),
          tags$tr(tags$td(strong("Compartment")), tags$td(row$compartment[1])),
          tags$tr(tags$td(strong("Formula")),     tags$td(row$formula[1])),
          tags$tr(tags$td(strong("Currency metabolite")), tags$td(if (row$is_currency[1]) "Yes" else "No"))
        )
      )
    } else {
      div(class = "text-muted", "Click a node in the network to see details here.")
    }
  })

  # REACTION TABLE ###################################

  output$gem_reaction_table <- renderDT({
    filtered_reactions() %>%
      arrange(distance, id) %>%
      select(id, name, distance, lower_bound, upper_bound, gene_reaction_rule, equation) %>%
      rename(
        `Reaction ID` = id,
        `Name`        = name,
        `Distance`    = distance,
        `Lower Bound` = lower_bound,
        `Upper Bound` = upper_bound,
        `GPR`         = gene_reaction_rule,
        `Equation`    = equation
      )
  }, rownames = FALSE, selection = "single", options = list(
    pageLength = 15, scrollX = TRUE
  ))

}

###################################################################################################
