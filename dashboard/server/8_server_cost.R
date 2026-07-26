###################################################################################################
# SERVER: Cost
###################################################################################################

page_cost_server <- function(input, output, session) {
  
  # DEFAULT DATA ------------------------------------------------------------------------------
  
  default_medium_df   <- read_csv("data/8_media.csv", show_col_types = FALSE)
  default_recovery_df <- read_csv("data/8_downstream.csv", show_col_types = FALSE)
  default_yield_df    <- read_csv("data/8_phb.csv", show_col_types = FALSE)
  
  # REACTIVE VALUES -----------------------------------------------------------------------------
  
  medium_data   <- reactiveVal(default_medium_df)
  recovery_data <- reactiveVal(default_recovery_df)
  yield_data    <- reactiveVal(default_yield_df)
  
  # MEDIUM COMPOSITION TABLE ---------------------------------------------------------------------
  
  output$cost_medium_table <- renderDT({
    datatable(
      medium_data(),
      rownames = FALSE,
      selection = "none",
      editable = list(target = "cell"),
      options = list(paging = FALSE, scrollX = TRUE, scrollY = "420px", dom = "t")
    )
  })
  
  medium_proxy <- dataTableProxy("cost_medium_table")
  
  observeEvent(input$cost_medium_table_cell_edit, {
    edit <- input$cost_medium_table_cell_edit
    df <- medium_data()
    
    col_name <- colnames(df)[edit$col + 1]
    new_val  <- edit$value
    
    if (col_name == "Component") {
      df[edit$row, edit$col + 1] <- new_val
    } else {
      suppressWarnings(num_val <- as.numeric(new_val))
      df[edit$row, edit$col + 1] <- ifelse(is.na(num_val), NA, num_val)
    }
    
    medium_data(df)
    replaceData(medium_proxy, df, resetPaging = FALSE, rownames = FALSE)
  })
  
  observeEvent(input$cost_add_component_row, {
    showModal(modalDialog(
      title = "Add New Component",
      textInput("cost_new_component_name", "Component Name", placeholder = "e.g. Glycerol"),
      numericInput("cost_new_component_cost", "Cost (PKR/kg)", value = 0, min = 0),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("cost_confirm_add_component", "Add", class = "btn-primary")
      )
    ))
  })
  
  observeEvent(input$cost_confirm_add_component, {
    req(input$cost_new_component_name)
    df <- medium_data()
    
    if (input$cost_new_component_name %in% df$Component) {
      showNotification("A component with that name already exists.", type = "error")
      return()
    }
    
    new_row <- as.list(rep(NA, ncol(df)))
    names(new_row) <- colnames(df)
    new_row[["Component"]] <- input$cost_new_component_name
    new_row[["Cost"]] <- ifelse(is.na(input$cost_new_component_cost), 0, input$cost_new_component_cost)
    
    df <- bind_rows(df, as_tibble(new_row))
    medium_data(df)
    removeModal()
  })
  
  observeEvent(input$cost_add_medium_col, {
    showModal(modalDialog(
      title = "Add New Medium",
      textInput("cost_new_medium_name", "Medium Name", placeholder = "e.g. RPB+Xylose"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("cost_confirm_add_medium", "Add", class = "btn-primary")
      )
    ))
  })
  
  observeEvent(input$cost_confirm_add_medium, {
    req(input$cost_new_medium_name)
    df <- medium_data()
    new_col_name <- input$cost_new_medium_name
    
    if (new_col_name %in% colnames(df)) {
      showNotification("A medium with that name already exists.", type = "error")
      return()
    }
    
    df[[new_col_name]] <- NA_real_
    medium_data(df)
    
    yd <- yield_data()
    if (!new_col_name %in% yd$`Medium Name`) {
      yd <- bind_rows(yd, tibble(`Medium Name` = new_col_name, `PHB` = 0.1))
      yield_data(yd)
    }
    
    removeModal()
  })
  
  # Keep medium choices synced
  observe({
    df <- medium_data()
    medium_choices <- setdiff(colnames(df), c("Component", "Cost"))
    current <- isolate(input$cost_select_medium)
    updateSelectInput(
      session, "cost_select_medium",
      choices = medium_choices,
      selected = if (!is.null(current) && current %in% medium_choices) current else medium_choices[1]
    )
  })
  
  # Keep plot media checkbox choices synced
  observe({
    df <- medium_data()
    medium_choices <- setdiff(colnames(df), c("Component", "Cost"))
    current <- isolate(input$cost_select_plot_media)
    
    selected_choices <- if (is.null(current)) medium_choices else intersect(current, medium_choices)
    
    updateCheckboxGroupInput(
      session, "cost_select_plot_media",
      choices = medium_choices,
      selected = selected_choices
    )
  })
  
  observeEvent(input$cost_plot_select_all, {
    df <- medium_data()
    medium_choices <- setdiff(colnames(df), c("Component", "Cost"))
    updateCheckboxGroupInput(session, "cost_select_plot_media", selected = medium_choices)
  })
  
  observeEvent(input$cost_plot_deselect_all, {
    updateCheckboxGroupInput(session, "cost_select_plot_media", selected = character(0))
  })
  
  # DOWNSTREAM RECOVERY TABLE --------------------------------------------------------------------
  
  output$cost_recovery_table <- renderDT({
    datatable(
      recovery_data(),
      rownames = FALSE,
      selection = "none",
      editable = list(target = "cell"),
      options = list(paging = FALSE, scrollX = TRUE, scrollY = "420px", dom = "t")
    )
  })
  
  recovery_proxy <- dataTableProxy("cost_recovery_table")
  
  observeEvent(input$cost_recovery_table_cell_edit, {
    edit <- input$cost_recovery_table_cell_edit
    df <- recovery_data()
    
    col_name <- colnames(df)[edit$col + 1]
    new_val  <- edit$value
    
    if (col_name %in% c("Method", "Description")) {
      df[edit$row, edit$col + 1] <- new_val
    } else {
      suppressWarnings(num_val <- as.numeric(new_val))
      df[edit$row, edit$col + 1] <- ifelse(is.na(num_val), NA, num_val)
    }
    
    recovery_data(df)
    replaceData(recovery_proxy, df, resetPaging = FALSE, rownames = FALSE)
  })
  
  observeEvent(input$cost_add_recovery_row, {
    showModal(modalDialog(
      title = "Add New Recovery Method",
      textInput("cost_new_recovery_method", "Method Name", placeholder = "e.g. NaOH Digestion"),
      textInput("cost_new_recovery_desc", "Description", placeholder = "Brief method description..."),
      numericInput("cost_new_recovery_cost", "Cost (PKR/1L)", value = 0, min = 0),
      numericInput("cost_new_recovery_scaling", "Scaling Factor", value = 1, min = 0, step = 0.1),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("cost_confirm_add_recovery", "Add Method", class = "btn-primary")
      )
    ))
  })
  
  observeEvent(input$cost_confirm_add_recovery, {
    req(input$cost_new_recovery_method)
    df <- recovery_data()
    
    if (input$cost_new_recovery_method %in% df$Method) {
      showNotification("A recovery method with that name already exists.", type = "error")
      return()
    }
    
    cost_col <- if ("Cost" %in% colnames(df)) "Cost" else "Cost (PKR/kg PHB)"
    desc_val <- if (is.null(input$cost_new_recovery_desc)) "" else input$cost_new_recovery_desc
    
    new_row <- as.list(rep(NA, ncol(df)))
    names(new_row) <- colnames(df)
    new_row[["Method"]] <- input$cost_new_recovery_method
    new_row[["Description"]] <- desc_val
    new_row[[cost_col]] <- ifelse(is.na(input$cost_new_recovery_cost), 0, input$cost_new_recovery_cost)
    new_row[["Scaling Factor"]] <- ifelse(is.na(input$cost_new_recovery_scaling), 1, input$cost_new_recovery_scaling)
    
    df <- bind_rows(df, as_tibble(new_row))
    recovery_data(df)
    removeModal()
  })
  
  # Keep recovery choices synced
  observe({
    df <- recovery_data()
    method_choices <- df$Method
    current <- isolate(input$cost_select_recovery)
    updateSelectInput(
      session, "cost_select_recovery",
      choices = method_choices,
      selected = if (!is.null(current) && current %in% method_choices) current else method_choices[1]
    )
  })
  
  # PHB YIELD TABLE -----------------------------------------------------------------------------
  
  output$cost_yield_table <- renderDT({
    datatable(
      yield_data(),
      rownames = FALSE,
      selection = "none",
      editable = list(target = "cell"),
      options = list(paging = FALSE, scrollX = TRUE, scrollY = "420px", dom = "t")
    )
  })
  
  yield_proxy <- dataTableProxy("cost_yield_table")
  
  observeEvent(input$cost_yield_table_cell_edit, {
    edit <- input$cost_yield_table_cell_edit
    df <- yield_data()
    
    col_name <- colnames(df)[edit$col + 1]
    new_val  <- edit$value
    
    if (col_name == "Medium Name") {
      df[edit$row, edit$col + 1] <- as.character(new_val)
    } else {
      suppressWarnings(num_val <- as.numeric(new_val))
      df[edit$row, edit$col + 1] <- ifelse(is.na(num_val), NA, num_val)
    }
    
    yield_data(df)
    replaceData(yield_proxy, df, resetPaging = FALSE, rownames = FALSE)
  })
  
  # IMPORT / EXPORT -------------------------------------------------------------------------------
  
  observeEvent(input$cost_import_media, {
    req(input$cost_import_media)
    imported <- tryCatch(
      read_csv(input$cost_import_media$datapath, show_col_types = FALSE),
      error = function(e) NULL
    )
    if (is.null(imported)) {
      showNotification("Failed to import Media CSV. Please check file format.", type = "error")
      return()
    }
    medium_data(imported)
    
    imp_media <- setdiff(colnames(imported), c("Component", "Cost"))
    yd <- yield_data()
    new_media <- setdiff(imp_media, yd$`Medium Name`)
    if (length(new_media) > 0) {
      yd <- bind_rows(yd, tibble(`Medium Name` = new_media, `PHB` = 0.1))
      yield_data(yd)
    }
    
    showNotification("Media composition imported successfully.", type = "message")
  })
  
  observeEvent(input$cost_import_recovery, {
    req(input$cost_import_recovery)
    imported <- tryCatch(
      read_csv(input$cost_import_recovery$datapath, show_col_types = FALSE),
      error = function(e) NULL
    )
    if (is.null(imported)) {
      showNotification("Failed to import Downstream Recovery CSV. Please check file format.", type = "error")
      return()
    }
    recovery_data(imported)
    showNotification("Downstream recovery methods imported successfully.", type = "message")
  })
  
  observeEvent(input$cost_import_yield, {
    req(input$cost_import_yield)
    imported <- tryCatch(
      read_csv(input$cost_import_yield$datapath, show_col_types = FALSE),
      error = function(e) NULL
    )
    if (is.null(imported)) {
      showNotification("Failed to import PHB Yield CSV. Please check file format.", type = "error")
      return()
    }
    yield_data(imported)
    showNotification("PHB yields imported successfully.", type = "message")
  })
  
  output$cost_export_media <- downloadHandler(
    filename = function() paste0("media_composition_", Sys.Date(), ".csv"),
    content = function(file) write_csv(medium_data(), file)
  )
  
  output$cost_export_recovery <- downloadHandler(
    filename = function() paste0("downstream_recovery_", Sys.Date(), ".csv"),
    content = function(file) write_csv(recovery_data(), file)
  )
  
  output$cost_export_yield <- downloadHandler(
    filename = function() paste0("phb_yield_", Sys.Date(), ".csv"),
    content = function(file) write_csv(yield_data(), file)
  )
  
  # COST CALCULATION & PLOT ------------------------------------------------------------------------
  
  medium_cost_per_liter <- reactive({
    req(input$cost_select_medium)
    df <- medium_data()
    req(input$cost_select_medium %in% colnames(df))
    
    amounts <- suppressWarnings(as.numeric(df[[input$cost_select_medium]]))
    costs   <- suppressWarnings(as.numeric(df[["Cost"]]))
    
    sum(amounts * costs / 1000, na.rm = TRUE)
  })
  
  phb_yield_val <- reactive({
    req(input$cost_select_medium)
    yd <- yield_data()
    row_match <- yd[yd$`Medium Name` == input$cost_select_medium, ]
    if (nrow(row_match) > 0) {
      col_name <- if ("PHB" %in% colnames(yd)) "PHB" else "PHB (g/1L volume)"
      val <- suppressWarnings(as.numeric(row_match[[col_name]][1]))
      if (!is.na(val) && val > 0) return(val)
    }
    return(1.0)
  })
  
  medium_cost_per_kg_phb <- reactive({
    (medium_cost_per_liter() / phb_yield_val()) * 1000
  })
  
  recovery_row <- reactive({
    req(input$cost_select_recovery)
    df <- recovery_data()
    req(input$cost_select_recovery %in% df$Method)
    df[df$Method == input$cost_select_recovery, ][1, ]
  })
  
  recovery_cost_per_kg_phb <- reactive({
    row <- recovery_row()
    req(row)
    col_name <- if ("Cost" %in% colnames(row)) "Cost" else "Cost (PKR/kg PHB)"
    cost_val <- suppressWarnings(as.numeric(row[[col_name]]))
    scaling  <- suppressWarnings(as.numeric(row[["Scaling Factor"]]))
    if (is.na(cost_val)) cost_val <- 0
    if (is.na(scaling))  scaling <- 1
    
    (cost_val / phb_yield_val()) * 1000 * scaling
  })
  
  total_cost_value <- reactive({
    medium_cost_per_kg_phb() + recovery_cost_per_kg_phb()
  })
  
  output$cost_medium_cost_out <- renderText({
    paste0(format(round(medium_cost_per_kg_phb(), 2), big.mark = ","), " PKR / kg PHB")
  })
  
  output$cost_recovery_cost_out <- renderText({
    paste0(format(round(recovery_cost_per_kg_phb(), 2), big.mark = ","), " PKR / kg PHB")
  })
  
  output$cost_total_cost_out <- renderText({
    paste0(format(round(total_cost_value(), 2), big.mark = ","), " PKR / kg PHB")
  })
  
  output$cost_comparison_plot <- renderPlotly({
    req(input$cost_select_recovery)
    
    selected_media <- input$cost_select_plot_media
    shiny::validate(
      need(!is.null(selected_media) && length(selected_media) > 0, "Please select at least one medium to display.")
    )
    
    m_df <- medium_data()
    y_df <- yield_data()
    r_df <- recovery_data()
    
    media_names <- intersect(selected_media, setdiff(colnames(m_df), c("Component", "Cost")))
    shiny::validate(
      need(length(media_names) > 0, "No valid media selected.")
    )
    
    rec_row <- r_df[r_df$Method == input$cost_select_recovery, ]
    cost_col <- if ("Cost" %in% colnames(r_df)) "Cost" else "Cost (PKR/kg PHB)"
    rec_cost_L <- if (nrow(rec_row) > 0) suppressWarnings(as.numeric(rec_row[[cost_col]][1])) else 0
    scaling    <- if (nrow(rec_row) > 0) suppressWarnings(as.numeric(rec_row[["Scaling Factor"]][1])) else 1
    if (is.na(rec_cost_L)) rec_cost_L <- 0
    if (is.na(scaling)) scaling <- 1
    
    comp_df <- tibble(
      Medium = media_names,
      Medium_Cost = sapply(media_names, function(m) {
        amounts <- suppressWarnings(as.numeric(m_df[[m]]))
        costs   <- suppressWarnings(as.numeric(m_df[["Cost"]]))
        m_cost_L <- sum(amounts * costs / 1000, na.rm = TRUE)
        
        y_row <- y_df[y_df$`Medium Name` == m, ]
        col_name <- if ("PHB" %in% colnames(y_df)) "PHB" else "PHB (g/1L volume)"
        y_val <- if (nrow(y_row) > 0) suppressWarnings(as.numeric(y_row[[col_name]][1])) else 1.0
        if (is.na(y_val) || y_val <= 0) y_val <- 1.0
        
        (m_cost_L / y_val) * 1000
      }),
      Recovery_Cost = sapply(media_names, function(m) {
        y_row <- y_df[y_df$`Medium Name` == m, ]
        col_name <- if ("PHB" %in% colnames(y_df)) "PHB" else "PHB (g/1L volume)"
        y_val <- if (nrow(y_row) > 0) suppressWarnings(as.numeric(y_row[[col_name]][1])) else 1.0
        if (is.na(y_val) || y_val <= 0) y_val <- 1.0
        
        (rec_cost_L / y_val) * 1000 * scaling
      })
    )
    
    plot_ly(comp_df, x = ~Medium, y = ~Medium_Cost, type = "bar", name = "Medium Cost") %>%
      add_trace(y = ~Recovery_Cost, name = "Downstream Recovery Cost") %>%
      layout(
        barmode = "stack",
        xaxis = list(title = "Medium"),
        yaxis = list(title = "Cost (PKR / kg PHB)"),
        legend = list(orientation = "h", x = 0.1, y = 1.1)
      )
  })
  
}

###################################################################################################