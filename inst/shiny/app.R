library(shiny)
library(bslib)
library(ggplot2)
library(plotly)
library(DT)
library(viralquantsim)

ui <- page_sidebar(
  title = "ViralQuantSim",
  theme = bs_theme(version = 5, bootswatch = "flatly"),

  sidebar = sidebar(
    h4("Assay configuration"),
    numericInput("n_samples", "Synthetic samples", 500, min = 50, max = 5000, step = 50),
    numericInput("seed", "Random seed", 42, min = 1, step = 1),
    sliderInput(
      "shared_recovery_sd",
      "Shared recovery variation (log10 SD)",
      min = 0,
      max = 0.8,
      value = 0.30,
      step = 0.05
    ),
    sliderInput(
      "ic_specific_recovery_sd",
      "IC-specific recovery noise (log10 SD)",
      min = 0,
      max = 0.8,
      value = 0.15,
      step = 0.05
    ),
    sliderInput(
      "qpcr_measurement_sd_ct",
      "qPCR measurement error (Ct SD)",
      min = 0.05,
      max = 1.0,
      value = 0.25,
      step = 0.05
    ),
    sliderInput(
      "target_library_specific_sd",
      "Library-specific variation (log10 SD)",
      min = 0,
      max = 1.0,
      value = 0.25,
      step = 0.05
    ),
    actionButton("simulate", "Run simulation", class = "btn-primary w-100"),
    hr(),
    downloadButton("download_data", "Download synthetic CSV", class = "w-100")
  ),

  navset_card_tab(
    nav_panel(
      "Overview",
      card(
        card_header("Scientific question"),
        p(
          "Under which technical conditions does a spike-in internal control improve ",
          "quantitative estimation in an integrated qPCR and sequencing workflow?"
        ),
        p(
          strong("Scope: "),
          "mechanistic simulation and methodological benchmarking. The application ",
          "does not contain participant-level data and is not a clinically validated calculator."
        )
      ),
      layout_columns(
        value_box("Samples", textOutput("sample_count")),
        value_box("qPCR detections", textOutput("detection_rate")),
        value_box("Sequencing dropouts", textOutput("dropout_rate"))
      )
    ),

    nav_panel(
      "Data",
      card(
        card_header("Synthetic observations"),
        DTOutput("data_table")
      )
    ),

    nav_panel(
      "qPCR",
      layout_columns(
        card(
          card_header("Target Ct versus true load"),
          plotlyOutput("target_ct_plot")
        ),
        card(
          card_header("Delta Ct versus true load"),
          plotlyOutput("delta_ct_plot")
        )
      )
    ),

    nav_panel(
      "Models",
      card(
        card_header("Model comparison"),
        DTOutput("model_table")
      ),
      card(
        card_header("Observed versus predicted"),
        selectInput(
          "selected_model",
          "Model",
          choices = c(
            "Target Ct" = "target_ct",
            "Delta Ct" = "delta_ct",
            "Reads" = "reads",
            "Reads + target Ct" = "reads_target_ct",
            "Reads + Delta Ct" = "reads_delta_ct"
          ),
          selected = "reads_delta_ct"
        ),
        plotlyOutput("prediction_plot")
      )
    ),

    nav_panel(
      "Methods",
      card(
        card_header("Model assumptions"),
        tags$ul(
          tags$li("The true target concentration is generated on a log10 scale."),
          tags$li("Target and internal control share some, but not all, technical variation."),
          tags$li("Ct is generated from an assay-specific intercept and slope, with censoring at the Ct limit."),
          tags$li("Sequencing counts use an overdispersed beta-binomial construction."),
          tags$li("Delta Ct is target Ct minus internal-control Ct."),
          tags$li("Model performance is reported both apparently and by repeated V-fold cross-validation.")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  simulation <- eventReactive(input$simulate, {
    parameters <- default_parameters()
    parameters$n_samples <- as.integer(input$n_samples)
    parameters$seed <- as.integer(input$seed)
    parameters$shared_recovery_sd <- input$shared_recovery_sd
    parameters$ic_specific_recovery_sd <- input$ic_specific_recovery_sd
    parameters$qpcr_measurement_sd_ct <- input$qpcr_measurement_sd_ct
    parameters$target_library_specific_sd <- input$target_library_specific_sd

    withProgress(message = "Simulating assay", value = 0.4, {
      result <- simulate_assay(parameters)
      incProgress(0.6)
      result
    })
  }, ignoreNULL = FALSE)

  candidate_models <- reactive({
    fit_candidate_models(simulation())
  })

  comparison <- reactive({
    compare_models(
      candidate_models(),
      simulation(),
      cross_validate = TRUE,
      v = 5,
      repeats = 3,
      seed = input$seed
    )
  })

  output$sample_count <- renderText({
    format(nrow(simulation()), big.mark = ",")
  })

  output$detection_rate <- renderText({
    paste0(round(mean(simulation()$target_detected) * 100, 1), "%")
  })

  output$dropout_rate <- renderText({
    paste0(round(mean(simulation()$sequencing_dropout) * 100, 1), "%")
  })

  output$data_table <- renderDT({
    visible <- simulation()[, !vapply(simulation(), is.list, logical(1)), drop = FALSE]
    datatable(
      visible,
      filter = "top",
      options = list(pageLength = 12, scrollX = TRUE)
    )
  })

  output$target_ct_plot <- renderPlotly({
    plot_data <- simulation()
    plot <- ggplot(plot_data, aes(true_log10_viral_load, target_ct)) +
      geom_point(alpha = 0.55) +
      geom_smooth(method = "lm", se = TRUE) +
      labs(x = "True viral load (log10 copies/mL)", y = "Target Ct") +
      theme_minimal()
    ggplotly(plot)
  })

  output$delta_ct_plot <- renderPlotly({
    plot_data <- simulation()
    plot <- ggplot(plot_data, aes(true_log10_viral_load, delta_ct)) +
      geom_point(alpha = 0.55) +
      geom_smooth(method = "lm", se = TRUE) +
      labs(x = "True viral load (log10 copies/mL)", y = "Delta Ct") +
      theme_minimal()
    ggplotly(plot)
  })

  output$model_table <- renderDT({
    result <- comparison()
    datatable(
      result,
      rownames = FALSE,
      options = list(pageLength = 5, scrollX = TRUE)
    )
  })

  output$prediction_plot <- renderPlotly({
    model <- candidate_models()[[input$selected_model]]
    ggplotly(plot_model_predictions(model, title = input$selected_model))
  })

  output$download_data <- downloadHandler(
    filename = function() {
      paste0("viralquantsim_", Sys.Date(), ".csv")
    },
    content = function(file) {
      visible <- simulation()[, !vapply(simulation(), is.list, logical(1)), drop = FALSE]
      utils::write.csv(visible, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
