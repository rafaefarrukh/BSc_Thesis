###################################################################################################
# SERVER: phaC
###################################################################################################

# Catalytic triad residue definitions (Bakta annotation + ColabFold geometric validation)
phac_triad_residues <- list(
  list(resi = 133, role = "Cys133 (nucleophile)",       color = "#e63946"),
  list(resi = 289, role = "Asp289 (orienting residue)", color = "#2a9d8f"),
  list(resi = 318, role = "His318 (general base)",      color = "#457b9d")
)

# Maps a representation-style choice to its m_style_*() constructor
phac_style_fn <- function(style_choice) {
  switch(style_choice,
    "Cartoon" = m_style_cartoon,
    "Stick"   = m_style_stick,
    "Sphere"  = m_style_sphere,
    "Line"    = m_style_line,
    m_style_cartoon
  )
}

# Maps a color-scheme choice to a style color value
phac_color_value <- function(color_choice) {
  if (identical(color_choice, "Spectrum (rainbow)")) "spectrum" else "#cfd8dc"
}

# Maps a background choice to its hex value
phac_bg_value <- function(bg_choice) {
  switch(bg_choice,
    "White" = "#ffffff",
    "Black" = "#000000",
    "Grey"  = "#333333",
    "#ffffff"
  )
}

# Applies triad highlight styling + labels to the live viewer given its output id.
# All m_*() calls here use the "id" (character) form of the API, which pushes a
# JS command straight to that existing widget instance -- there is no proxy object.
phac_highlight_triad <- function(id) {
  for (res in phac_triad_residues) {
    m_add_style(
      id    = id,
      sel   = m_sel(resi = res$resi),
      style = m_style_stick(color = res$color, radius = 0.3)
    )
    m_add_res_labels(
      id    = id,
      sel   = m_sel(resi = res$resi),
      style = m_style_label(
        backgroundColor = res$color,
        fontColor = "white",
        fontSize = 12,
        showBackground = TRUE
      )
    )
  }
}

# Applies the current representation + color scheme selections to the whole model
phac_apply_base_style <- function(id, style_choice, color_choice) {
  style_fn <- phac_style_fn(style_choice)
  m_set_style(id = id, style = style_fn(color = phac_color_value(color_choice)))
}

page_phac_server <- function(input, output, session, pdb_path = "data/phaC.pdb") {

  pdb_text <- paste(readLines(pdb_path, warn = FALSE), collapse = "\n")

  output$phac_structure_viewer <- renderR3dmol({
    r3dmol() |>
      m_add_model(data = pdb_text, format = "pdb") |>
      m_set_style(style = m_style_cartoon(color = "#cfd8dc")) |>
      m_zoom_to()
  })

  # --- Representation / color scheme ---------------------------------------
  observeEvent(list(input$phac_style_select, input$phac_color_select), {
    req(input$phac_style_select, input$phac_color_select)
    m_remove_all_labels(id = "phac_structure_viewer")
    phac_apply_base_style(
      id = "phac_structure_viewer",
      style_choice = input$phac_style_select,
      color_choice = input$phac_color_select
    )
  }, ignoreInit = TRUE)

  # --- Background ------------------------------------------------------------
  observeEvent(input$phac_bg_select, {
    m_set_background_color(id = "phac_structure_viewer", hex = phac_bg_value(input$phac_bg_select))
  }, ignoreInit = TRUE)

  # --- Spin toggle -------------------------------------------------------------
  observeEvent(input$phac_spin_toggle, {
    if (isTRUE(input$phac_spin_toggle)) {
      m_spin(id = "phac_structure_viewer", axis = "y", speed = 1)
    } else {
      m_spin(id = "phac_structure_viewer", axis = "y", speed = 0)
    }
  }, ignoreInit = TRUE)

  # --- Zoom in / out -----------------------------------------------------------
  observeEvent(input$phac_btn_zoom_in, {
    m_zoom(id = "phac_structure_viewer", factor = 1.3, animationDuration = 300)
  })

  observeEvent(input$phac_btn_zoom_out, {
    m_zoom(id = "phac_structure_viewer", factor = 0.7, animationDuration = 300)
  })

  # --- Reset view --------------------------------------------------------------
  observeEvent(input$phac_btn_reset, {
    updateSelectInput(session, "phac_style_select", selected = "Cartoon")
    updateSelectInput(session, "phac_color_select", selected = "Uniform (grey)")
    updateSelectInput(session, "phac_bg_select", selected = "White")
    updateCheckboxInput(session, "phac_spin_toggle", value = TRUE)

    m_remove_all_labels(id = "phac_structure_viewer")
    m_set_style(id = "phac_structure_viewer", style = m_style_cartoon(color = "#cfd8dc"))
    m_set_background_color(id = "phac_structure_viewer", hex = "#ffffff")
    m_spin(id = "phac_structure_viewer", axis = "y", speed = 0)
    m_zoom_to(id = "phac_structure_viewer")
  })

  # --- Catalytic triad exploration ----------------------------------------------
  observeEvent(input$phac_btn_full, {
    m_remove_all_labels(id = "phac_structure_viewer")
    phac_apply_base_style(
      id = "phac_structure_viewer",
      style_choice = input$phac_style_select,
      color_choice = input$phac_color_select
    )
  })

  observeEvent(input$phac_btn_triad, {
    m_remove_all_labels(id = "phac_structure_viewer")
    phac_apply_base_style(
      id = "phac_structure_viewer",
      style_choice = input$phac_style_select,
      color_choice = input$phac_color_select
    )
    phac_highlight_triad(id = "phac_structure_viewer")
  })

  observeEvent(input$phac_btn_zoom, {
    m_remove_all_labels(id = "phac_structure_viewer")
    phac_apply_base_style(
      id = "phac_structure_viewer",
      style_choice = input$phac_style_select,
      color_choice = input$phac_color_select
    )
    phac_highlight_triad(id = "phac_structure_viewer")
    m_zoom_to(id = "phac_structure_viewer", sel = m_sel(resi = c(133, 289, 318)))
  })
}

###################################################################################################
