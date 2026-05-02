library(shiny)
library(igraph)
library(Matrix)

for (r_file in sort(list.files("PandemicsR-main/R", pattern = "\\.R$", full.names = TRUE))) {
  sys.source(r_file, envir = environment())
}

max_supported_n <- 200000L
max_supported_m <- 1000L
max_supported_t <- 500L
exact_simulation_node_limit <- 500L
graph_node_limit <- 250L

resolve_app_password <- function() {
  password <- Sys.getenv("APP_PASSWORD", unset = "")
  if (nzchar(password)) {
    return(password)
  }

  "voter"
}

password_gate_card <- function(error_message = NULL) {
  div(
    class = "password-gate-wrap",
    div(
      class = "password-gate-card",
      tags$h3("Protected Demo"),
      tags$p("Enter the shared password to open the dashboard."),
      passwordInput("accessPassword", "Password"),
      actionButton("unlockApp", "Open Dashboard"),
      if (!is.null(error_message) && nzchar(error_message)) {
        tags$p(class = "password-gate-error", error_message)
      }
    )
  )
}

param_help <- function(text) {
  tags$p(class = "param-help", text)
}

#==========================
# --- UI ---
#==========================
dashboard_ui <- tagList(
  titlePanel(
    "Voter - Schelling - Epidemic Multi-Membership Simulation",
    windowTitle = "PandemicsR Sandbox"
  ),

  sidebarLayout(

    sidebarPanel(
      numericInput("n", "Number of individuals", value = 3000, min = 5, max = max_supported_n),
      param_help("Total number of individuals in the simulation."),
      numericInput("m", "Number of groups", value = 3, min = 2, max = max_supported_m),
      param_help("Total number of groups that individuals can join."),
      numericInput("timesteps", "Iterations", value = 100, min = 1, max = max_supported_t),
      param_help("Maximum Gillespie time horizon for one simulation run."),

      sliderInput("lambda", "RIG weight parameter lambda", min = 0, step = 0.01, max = 100, value = 0.5),
      param_help("Higher values create a denser initial individual-group membership structure."),
      sliderInput("c_param", "Schelling: Edge addition rate param c", min = 0, step = 0.01, max = 1, value = 0.4),
      param_help("Controls how quickly outsiders join groups."),

      sliderInput("gamma", "Voter: light-state gamma", min = 0, step = 0.01, max = 20, value = 0.1),
      param_help("Baseline opinion-change rate for the light shades."),
      sliderInput("gamma_dark", "Voter: dark-state gamma", min = 0, step = 0.01, max = 5, value = 0.02),
      param_help("Dark shades are more reluctant to change opinion."),
      sliderInput("infected_dark_multiplier", "Voter: infected dark-state multiplier", min = 1, step = 0.1, max = 10, value = 1.5),
      param_help("Boosts dark-state opinion changes when those individuals are infected."),

      sliderInput("beta_plus", "Schelling: beta_plus", min = 0, step = 0.01, max = 1, value = 0.5),
      param_help("Leaving rate when a camp is underrepresented inside a group."),
      sliderInput("beta_minus", "Schelling: beta_minus", min = 0, step = 0.01, max = 1, value = 0.2),
      param_help("Leaving rate when a camp is already well represented inside a group."),
      sliderInput("T_threshold", "Schelling: T_threshold", min = 0, step = 0.01, max = 1, value = 0.3),
      param_help("Minimum same-camp share needed to avoid the higher leaving rate."),

      sliderInput("Numopinions", "Number of Opinions", min = 2, step = 2, max = 4, value = 4),
      param_help("Choose either the legacy 2-state mode or the ordered 4-state mode."),
      sliderInput("beta_red_epi", "Epidemic: infection rate red camp", min = 0, step = 0.01, max = 2, value = 1.1),
      param_help("Dark red and light red share this infection rate in the random-mixing SIR process."),
      sliderInput("beta_blue_epi", "Epidemic: infection rate blue camp", min = 0, step = 0.01, max = 2, value = 0.18),
      param_help("Dark blue and light blue share this infection rate."),
      sliderInput("gamma_sir", "Epidemic: recovery rate", min = 0, step = 0.01, max = 2, value = 0.30),
      param_help("Common recovery rate for infected individuals."),
      sliderInput("initial_infected_fraction", "Epidemic: initial infected fraction", min = 0, step = 0.01, max = 0.5, value = 0.05),
      param_help("Fraction of individuals seeded as infected at time 0."),

      checkboxInput("runVoter", "Enable voter dynamics", value = TRUE),
      param_help("Turns opinion updates on or off."),
      checkboxInput("runSchelling", "Enable Schelling dynamics", value = TRUE),
      param_help("Turns group joining and leaving dynamics on or off."),
      actionButton("runSim", "Run Simulation"),
      uiOutput("runStatus"),
      checkboxInput("show_rig0", "Show initial Graph", value = FALSE),
      param_help("Displays the initial RIG and bipartite graph for small runs."),
      checkboxInput("show_rig", "Show final Graph", value = FALSE),
      param_help("Displays the final RIG and bipartite graph for small runs.")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Graphs",
          plotOutput("rig0Plot", height = "500px"),
          plotOutput("rigPlot", height = "500px"),
          plotOutput("bipartite0Plot", height = "500px"),
          plotOutput("bipartitePlot", height = "500px")
        ),
        tabPanel("Opinion dynamics",
          h3("Overall opinions"),
          plotOutput("histo"),

          h3("Time evolution of opinions"),
          plotOutput("fracPlot", width = "650px", height = "400px")
        ),
        tabPanel("Epidemics",
          h3("SIR dynamics"),
          plotOutput("sirPlot", width = "650px", height = "400px"),
          h3("SIR dynamics by camp"),
          plotOutput("campSirTimePlot", width = "650px", height = "520px"),
          h3("Cumulative infections by camp"),
          plotOutput("infectionPlot", width = "650px", height = "360px"),
          h3("Final attack rate by camp"),
          plotOutput("campAttackPlot", width = "500px", height = "320px"),
          h3("Final SIR counts by camp"),
          tableOutput("campSirTable"),
          verbatimTextOutput("epiSummary")
        ),
        tabPanel("Others",
          h3("Group-wise initial opinions"),
          plotOutput("histogramGroup0"),

          h3("Group-wise final opinions"),
          plotOutput("histogramGroup")
        )
      )
    )
  )
)

ui <- fluidPage(
  title = "PandemicsR Sandbox",
  tags$head(
    tags$style(HTML("
      .password-gate-wrap {
        display: flex;
        justify-content: center;
        margin-top: 32px;
      }
      .password-gate-card {
        width: 100%;
        max-width: 420px;
        padding: 24px;
        background: #ffffff;
        border: 1px solid #d0d7de;
        border-radius: 16px;
        box-shadow: 0 12px 36px rgba(15, 23, 42, 0.08);
      }
      .password-gate-error {
        margin-top: 12px;
        color: #b42318;
        font-weight: 600;
      }
      .param-help {
        margin-top: -8px;
        margin-bottom: 12px;
        color: #5b6879;
        font-size: 12px;
        line-height: 1.35;
      }
      .run-status {
        margin-top: 8px;
        margin-bottom: 12px;
        font-size: 12px;
        line-height: 1.35;
      }
      .run-status-pending {
        color: #9b1c1c;
      }
      .run-status-ok {
        color: #1d6f42;
      }
      .run-status-info {
        color: #5b6879;
      }
    "))
  ),
  uiOutput("appShell")
)

#==========================
# --- Server ---
#==========================
server <- function(input, output, session) {
  is_authenticated <- reactiveVal(FALSE)
  auth_error <- reactiveVal(NULL)
  sim_result <- reactiveVal(NULL)
  sim_error <- reactiveVal(NULL)
  initial_sim_done <- reactiveVal(FALSE)
  last_run_signature <- reactiveVal(NULL)
  expected_password <- resolve_app_password()
  required_input_ids <- c(
    "n", "m", "timesteps", "lambda", "c_param", "gamma",
    "gamma_dark", "infected_dark_multiplier",
    "beta_plus", "beta_minus", "T_threshold", "Numopinions",
    "beta_red_epi", "beta_blue_epi", "gamma_sir", "initial_infected_fraction",
    "runVoter", "runSchelling"
  )

  inputs_ready <- function() {
    all(vapply(required_input_ids, function(id) !is.null(input[[id]]), logical(1)))
  }

  current_signature <- reactive({
    req(inputs_ready())
    vapply(required_input_ids, function(id) paste(input[[id]], collapse = ","), character(1))
  })

  validate_scalability_inputs <- function() {
    n_value <- as.integer(input$n)
    m_value <- as.integer(input$m)
    t_value <- as.numeric(input$timesteps)

    if (!is.finite(n_value) || n_value < 5L || n_value > max_supported_n) {
      stop(sprintf("Number of individuals must be between 5 and %s.", max_supported_n), call. = FALSE)
    }
    if (!is.finite(m_value) || m_value < 2L || m_value > max_supported_m || m_value > n_value) {
      stop(sprintf("Number of groups must be between 2 and %s, and no larger than the number of individuals.", max_supported_m), call. = FALSE)
    }
    if (!is.finite(t_value) || t_value < 1 || t_value > max_supported_t) {
      stop(sprintf("Iterations must be between 1 and %s.", max_supported_t), call. = FALSE)
    }

    invisible(TRUE)
  }

  output$appShell <- renderUI({
    if (isTRUE(is_authenticated())) {
      dashboard_ui
    } else {
      password_gate_card(auth_error())
    }
  })

  output$runStatus <- renderUI({
    req(isTRUE(is_authenticated()))
    if (!inputs_ready()) {
      return(NULL)
    }

    last_signature <- last_run_signature()
    if (is.null(last_signature)) {
      return(tags$p(class = "run-status run-status-info", "The dashboard runs once when it opens. Use Run Simulation after changing parameters."))
    }

    if (!identical(last_signature, current_signature())) {
      return(tags$p(class = "run-status run-status-pending", "Parameters changed. Press Run Simulation to refresh all plots."))
    }

    tags$p(class = "run-status run-status-ok", "Plots match the current settings.")
  })

  observeEvent(input$unlockApp, {
    entered_password <- input$accessPassword
    if (is.null(entered_password)) {
      entered_password <- ""
    }
    if (identical(entered_password, expected_password)) {
      auth_error(NULL)
      sim_error(NULL)
      sim_result(NULL)
      initial_sim_done(FALSE)
      last_run_signature(NULL)
      is_authenticated(TRUE)
    } else {
      auth_error("Incorrect password.")
    }
  }, ignoreInit = TRUE)

  observe({
    req(inputs_ready())
    if (input$n > graph_node_limit) {
      updateCheckboxInput(session, "show_rig0", value = FALSE)
      updateCheckboxInput(session, "show_rig", value = FALSE)
    }
  })

  run_current_simulation <- function() {
    if (!inputs_ready()) {
      return(invisible(FALSE))
    }

    result <- tryCatch(
      withProgress(message = "Running simulation", value = 0.15, {
        validate_scalability_inputs()
        simulate_hybrid_model(
          n = input$n,
          m = input$m,
          t_max = input$timesteps,
          lambda = input$lambda,
          c_param = input$c_param,
          gamma_light = input$gamma,
          gamma_dark = input$gamma_dark,
          infected_dark_multiplier = input$infected_dark_multiplier,
          beta_plus = input$beta_plus,
          beta_minus = input$beta_minus,
          T_threshold = input$T_threshold,
          num_opinions = input$Numopinions,
          run_voter = isTRUE(input$runVoter),
          run_schelling = isTRUE(input$runSchelling),
          beta_red = input$beta_red_epi,
          beta_blue = input$beta_blue_epi,
          gamma_sir = input$gamma_sir,
          initial_infected_fraction = input$initial_infected_fraction,
          simulation_mode = "auto",
          exact_threshold = exact_simulation_node_limit,
          graph_threshold = graph_node_limit
        )
      }),
      error = function(err) err
    )

    if (inherits(result, "error")) {
      sim_error(conditionMessage(result))
      sim_result(NULL)
      return(invisible(FALSE))
    }

    sim_error(NULL)
    sim_result(result)
    last_run_signature(isolate(current_signature()))
    invisible(TRUE)
  }

  observe({
    req(isTRUE(is_authenticated()))
    req(!isTRUE(initial_sim_done()))
    req(inputs_ready())
    if (isTRUE(run_current_simulation())) {
      initial_sim_done(TRUE)
    }
  })

  observeEvent(input$runSim, {
    req(isTRUE(is_authenticated()))
    if (isTRUE(run_current_simulation())) {
      initial_sim_done(TRUE)
    }
  }, ignoreInit = TRUE)

  simData <- reactive({
    req(isTRUE(is_authenticated()))
    if (!is.null(sim_error())) {
      validate(need(FALSE, sim_error()))
    }
    req(sim_result())
    sim_result()
  })

  #==========================
  # --- Output ---
  #==========================
  
  output$rig0Plot <- renderPlot({
    req(isTRUE(is_authenticated()))
    req(simData())
    req(input$show_rig0, cancelOutput = TRUE)
    req(isTRUE(simData()$graph_available), cancelOutput = TRUE)
    
    visual_step_multi(
      simData()$RIG0,
      simData()$opinion_history[,1],
      simData()$num_opinions
    )
  })
  
  output$rigPlot <- renderPlot({
    req(isTRUE(is_authenticated()))
    req(simData())
    req(input$show_rig, cancelOutput = TRUE)
    req(isTRUE(simData()$graph_available), cancelOutput = TRUE)
    
    visual_step_multi(
      simData()$RIG,
      simData()$opinion_history[, ncol( simData()$opinion_history)],
      simData()$num_opinions
    )
  })
  
  output$bipartite0Plot <- renderPlot({
    req(isTRUE(is_authenticated()))
    req(simData())
    req(input$show_rig0, cancelOutput = TRUE)
    req(isTRUE(simData()$graph_available), cancelOutput = TRUE)

    visual_bipartite(simData()$B0, simData()$opinion_history[,1], simData()$num_opinions)
  })
  
  output$bipartitePlot <- renderPlot({
    req(isTRUE(is_authenticated()))
    req(simData())
    req(input$show_rig, cancelOutput = TRUE)
    req(isTRUE(simData()$graph_available), cancelOutput = TRUE)

    visual_bipartite(simData()$B, simData()$opinion_history[,ncol( simData()$opinion_history)], simData()$num_opinions)
  })

  output$fracPlot <- renderPlot({ req(isTRUE(is_authenticated())); req(simData()); visual_step_time(simData()$frac_mat, simData()$num_opinions, simData()$time_history) })
  output$histo <- renderPlot({
    req(isTRUE(is_authenticated()))
    req(simData())
    if (identical(simData()$simulation_mode, "aggregate")) {
      visual_opinion_shares(simData()$frac_mat, simData()$num_opinions)
    } else {
      visual_histo(simData()$opinion_history, simData()$num_opinions)
    }
  })
  output$sirPlot <- renderPlot({ req(isTRUE(is_authenticated())); req(simData()); visual_sir_time(simData()$sir_mat, simData()$time_history) })
  output$campSirTimePlot <- renderPlot({ req(isTRUE(is_authenticated())); req(simData()); visual_sir_camp_time(simData()$camp_sir_history, simData()$time_history) })
  output$infectionPlot <- renderPlot({ req(isTRUE(is_authenticated())); req(simData()); visual_cumulative_infections(simData()$infection_events, simData()$final_time, simData()$camp_labels) })
  output$campAttackPlot <- renderPlot({ req(isTRUE(is_authenticated())); req(simData()); visual_epidemic_camps(simData()$final_camp_sir) })
  output$campSirTable <- renderTable({
    req(isTRUE(is_authenticated()))
    req(simData())
    as.data.frame(round(simData()$final_camp_sir, if (identical(simData()$simulation_mode, "aggregate")) 1L else 0L))
  }, rownames = TRUE)
  output$epiSummary <- renderPrint({
    req(isTRUE(is_authenticated()))
    req(simData())
    cat(sprintf("Simulation mode: %s\n", simData()$simulation_mode))
    cat(sprintf("Mode note: %s\n", simData()$model_notes))
    cat(sprintf("Overall attack rate: %.2f\n", simData()$overall_attack_rate))
    cat(sprintf("Red camp attack rate: %.2f\n", simData()$camp_attack_rate[["Red camp"]]))
    cat(sprintf("Blue camp attack rate: %.2f\n", simData()$camp_attack_rate[["Blue camp"]]))
    cat(sprintf("Final recorded simulation time: %.2f\n", simData()$final_time))
  })
  output$histogramGroup0 <- renderPlot({
    req(isTRUE(is_authenticated()))
    req(simData())
    req(!identical(simData()$simulation_mode, "aggregate"), cancelOutput = TRUE)
    visual_histo_pergroup(simData()$opinion_history[,1], simData()$num_opinions, simData()$members0)
  })
  output$histogramGroup <- renderPlot({
    req(isTRUE(is_authenticated()))
    req(simData())
    req(!identical(simData()$simulation_mode, "aggregate"), cancelOutput = TRUE)
    visual_histo_pergroup(simData()$opinion_history[,ncol( simData()$opinion_history)], simData()$num_opinions, simData()$members)
  })
  
}

#==========================
# --- Run App ---
#==========================
app <- shinyApp(ui, server)
app
