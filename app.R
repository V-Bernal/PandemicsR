library(shiny)
library(igraph)
library(Matrix)
library(PandemicsR)
library(zip)

#==========================
# --- UI ---
#==========================
ui <- fluidPage(

  #==========================
  # Title
  #==========================
  titlePanel(h3("Voter - Schelling Multi-Membership with Epidemics")),


  #==========================
  # Side panel
  #==========================
  sidebarLayout(

    sidebarPanel(

      # Section 1: Parameters
      wellPanel(
        h4("Network parameters"),
        fluidRow(
          column(6, numericInput( "n", "Individuals", value = 15) ),
          column(6, numericInput("m", "Groups", value = 3) )
          ),
        
        fluidRow(column(6, numericInput("timesteps", "Max time", value = 100)),
          column(6,numericInput("lambda", "RIG weight (lambda) ",value = 1)) )
              ),
      
      # wellPanel(
      #   radioButtons(
      #     "parameter_mode",
      #     "Simulation scenarios:",
      #     choices = c(
      #       "Manual" = "manual",
      #       "Predefined scenario" = "scenario"
      #     ),
      #     selected = "manual"
      #   ),
      #   
      #   conditionalPanel(
      #     condition = "input.parameter_mode == 'scenario'",
      #   radioButtons(
      #   "scenario",
      #   "Simulation scenario:",
      #   choices = c(
      #     "Custom" = "custom",
      #     "1. Resilient / Low-risk" = "resilient",
      #     "2. Polarized / Segregated" = "polarized",
      #     "3. Radicalization-dominated" = "radicalization",
      #     "4. Epidemic-dominated" = "epidemic"
      #   ),
      #   selected = "resilient"
      # ),
      # textOutput("scenario_description")
      #   ),
      # 
      # 
      # ),
      # Section 2: Schelling
      wellPanel(
        h4("Schelling model"),
        checkboxInput("runSchelling", "Schelling model", value = FALSE),
        
        conditionalPanel(
          condition = "input.runSchelling",
          
        sliderInput("c_param", "join group rate c", min = 0, step = 0.1, max = 1000, value = 1),
        sliderInput("beta_plus", "leave rate above tolerance", min = 0, step = 0.01, max = 1, value = 0.5),
        sliderInput("beta_minus", "leave rate below tolerance", min = 0, step = 0.01, max = 1, value = 0.2),
        sliderInput("T_threshold", "tolerance threshold", min = 0, step = 0.01, max = 1, value = 0.3)        )
        ),

      # Section 3: Voter
      wellPanel(
        h4("Voter's model"),
        checkboxInput("runVoter", "Voter model", value = FALSE),
        sliderInput(
          "Numopinions",
          "Number of Opinions",
          min = 2,
          step = 2,
          max = 4,
          value = 4
        ),
        conditionalPanel(
          condition = "input.runVoter",
        sliderInput("gamma", "Opinions rate", min = 0, step = 0.1, max = 100, value = 5),
        #sliderInput("Numopinions", "Number of Opinions", min = 2, step = 2, max = 4, value = 4),

        # Section 4: Extremes
        #h5("Extremes"),
        sliderInput("alpha", "radicalization rate", min = 0, step = 0.1, max = 1000, value = 1),
        sliderInput("alpha_deradicalization", "deradicalization rate*", min = 0, step = 0.1, max = 1000, value = 1),
        h6("*stubborn opinions: set de-radicalization to zero"),
        sliderInput("alpha0_rad", "spontaneous radicalization rate", min = 0, step = 0.1, max = 1000, value = 0),
        sliderInput("alpha0_derad", "spontaneous deradicalization rate", min = 0, step = 0.1, max = 1000, value = 0)
        )
        ),

      # Section 5: Epidemics
      wellPanel(
        h4("Epidemics model"),
        checkboxInput("runEpidemic", "Epidemic model", value = FALSE),
        
        conditionalPanel(
          
        condition = "input.runEpidemic",
        sliderInput("I0", "Fraction of Infected", min = 0.01, step = 0.1, max = 1, value = 0.1),
        sliderInput("gamma_epi", "Recovery rate", min = 0.00, step = 0.1, max = 1000, value = 1),
        sliderInput("beta_red_red", "Infection rate red-red", min = 0, step = 0.1, max = 1000, value = 1),
        sliderInput("beta_blue_blue", "Infection rate blue-blue", min = 0, step = 0.1, max = 1000, value = 1),
        sliderInput("beta_red_blue", "Infection rate red-blue", min = 0, step = 0.1, max = 1000, value = 1),
        sliderInput("beta_blue_red", "Infection rate blue-red", min = 0, step = 0.1, max = 1000, value = 1)
        
        ),
      
      # Section 5.1: Epidemics activation
      wellPanel(
        h4("Epidemic activation"),
        
        selectInput(
          "epi_trigger",
          "Start epidemic after:",
          choices = c(
            "Fixed time" = "time",
            "Opinion and Segregation stability" = "global_stability"
          ),
          selected = "time"
        ),
        
        conditionalPanel(
          condition = "input.epi_trigger == 'time'",
          
          numericInput(
            "epi_time",
            "Epidemic start time:",
            value = 0,
            min = 0,
            step = 1
          )
        ),
        
        conditionalPanel(
          condition = "input.epi_trigger == 'global_stability'",
          
          numericInput(
            "stability_window",
            "Stability window (time):",
            value = 10,
            min = 1,
            step = 1
          ),
          
          numericInput(
            "stability_threshold",
            "Stability threshold:",
            value = 0.01,
            min = 0,
            max = 1,
            step = 0.01
          )
        )
      )

    )
      ,
      # Section 6: Visualization
      wellPanel(
        h4("Visualization"),
        checkboxInput("show_rig0", "Show initial Graph", value = FALSE),
        checkboxInput("show_rig", "Show final Graph", value = FALSE)
      ),

      # Run
      actionButton("runSim", "Run Simulation"),
      downloadButton("downloadAllPlots", "Download All Plots")
    ),

    #==========================
    # Main panel: Visualizations
    #==========================
    mainPanel(

      tabsetPanel(
        tabPanel("Graphs",

                 h4("Stopping condition"),
                 textOutput("stopReason"),
                 textOutput("compTime"),
                 textOutput("numEvents"),

                 fluidRow(
                   column(6,
                          h4("Initial RIG"),
                          plotOutput("rig0Plot", height = "400px")
                   ),
                   column(6,
                          h4("Final RIG"),
                          plotOutput("rigPlot", height = "400px")
                   )
                 ),

                 h4("Bipartite graph"),
                 p("Connections between individuals and groups."),

                 fluidRow(
                   column(6,
                          h4("Initial bipartite"),
                          plotOutput("bipartite0Plot", height = "400px")
                   ),
                   column(6,
                          h4("Final bipartite"),
                          plotOutput("bipartitePlot", height = "400px")
                   )
                 )
        ),

        tabPanel("Voter's dynamics",
                 h4("Time evolution of opinions"),
                 p("Fraction of individuals holding each opinion state over time."),
                 plotOutput("fracPlot", width = "500px", height = "400px"),

                 h4("Overall opinions"),
                 p("Distribution of opinions across the entire population."),
                 plotOutput("histo"),

                 h4("Group-wise initial opinions"),
                 p("Opinion distribution within each group at the start."),
                 plotOutput("histogramGroup0"),

                 h4("Group-wise final opinions"),
                 p("Opinion distribution within each group at the end of the simulation."),
                 plotOutput("histogramGroup")



        ),
        
        tabPanel("Epidemic's dynamics",
                 
                 h4("Starting condition epidemics"),
                 textOutput("startReasonEpi"),


                 h4("SIR overall"),
                 p("Total population fractions in Susceptible (S), Infected (I), and Recovered (R) over time."),
                 plotOutput("SIR1"),

                 h4("SIR red opinion"),
                 p("SIR dynamics restricted to individuals currently in the red opinion camp."),
                 plotOutput("SIR2"),

                 h4("SIR blue opinion"),
                 p("SIR dynamics restricted to individuals currently in the blue opinion camp."),
                 plotOutput("SIR3"),

                 h4("SIR members"),
                 p("SIR dynamics for individuals belonging to at least one group."),
                 plotOutput("SIR4"),

                 h4("SIR isolated"),
                 p("SIR dynamics for individuals not belonging to any group."),
                 plotOutput("SIR5")
        )
      )
    )
  )
)

#==========================
# --- Server ---
#==========================
server <- function(input, output, session) {

  simData <- eventReactive(input$runSim, {
    
    p <- list(
      
      # ==========================
      # Network
      # ==========================
      
      n = input$n,
      m = input$m,
      t_max = input$timesteps,
      lambda = input$lambda,
      
      # ==========================
      # Model switches
      # ==========================
      
      runVoter = isTRUE(input$runVoter),
      runSchelling = isTRUE(input$runSchelling),
      runEpidemic = isTRUE(input$runEpidemic),
      
      # ==========================
      # Opinions
      # ==========================
      
      num_opinions = input$Numopinions,
      
      # ==========================
      # Schelling
      # ==========================
      
      c_param = input$c_param,
      beta_plus = input$beta_plus,
      beta_minus = input$beta_minus,
      T_threshold = input$T_threshold,
      
      # ==========================
      # Voter / radicalization
      # ==========================
      
      gamma = input$gamma,
      
      alpha = input$alpha,
      alpha_deradicalization =
        input$alpha_deradicalization,
      
      alpha0_rad = input$alpha0_rad,
      alpha0_derad = input$alpha0_derad,
      
      # ==========================
      # Epidemic
      # ==========================
      
      beta_red_red = input$beta_red_red,
      beta_red_blue = input$beta_red_blue,
      beta_blue_red = input$beta_blue_red,
      beta_blue_blue = input$beta_blue_blue,
      
      gamma_epi = input$gamma_epi,
      I0 = input$I0,
      
      # ==========================
      # Epidemic activation
      # ==========================
      
      epi_trigger = input$epi_trigger,
      
      epi_time =
        if (identical(input$epi_trigger, "time"))
          input$epi_time
      else
        NULL,
      
      stability_window =
        if (identical(input$epi_trigger, "global_stability"))
          input$stability_window
      else
        NULL,
      
      stability_threshold =
        if (identical(input$epi_trigger, "global_stability"))
          input$stability_threshold
      else
        NULL
    )
    
    # print("PARAMETERS:")
    # print(p)
    run_simulation(p)
  })

  #==========================
  # Section 6: Visualization
  #==========================
  
  # Stopping reason
  output$stopReason <- renderText({
    req(simData())
    paste("Simulation stopped:", simData()$stop_reason,
          "| Final time:", round(simData()$final_time, 2))
  })

  output$compTime <- renderText({
    req(simData())
    paste0("Computation time: ",
           round(as.numeric(simData()$comp_time, units = "secs"), 3),
           " seconds")
  })

  output$numEvents <- renderText({
    req(simData())
    paste0("Number of events: ",
           round(as.numeric(simData()$event_counter, units = "secs"), 3),
           " events")
  })

  output$scenario_description <- renderText({
    
    switch(
      input$scenario,
      
      resilient =
        "Moderate social dynamics, lower transmission and faster recovery.",
      
      polarized =
        "Strong within-camp dynamics and segregation; cross-camp transmission is lower.",
      
      radicalization =
        "Radicalization dominates deradicalization and epidemic persistence is favored.",
      
      epidemic =
        "Strong epidemic pressure with active social feedback."
    )
  })
  # Network 
  output$rig0Plot <- renderPlot({
    req(simData())           
    req(input$show_rig0)     

    visual_step_multi(
      simData()$RIG0,
      simData()$opinion_history[,1],
      simData()$num_opinions
    )
  })

  output$rigPlot <- renderPlot({
    req(simData())           
    req(input$show_rig)     

    visual_step_multi(
      simData()$RIG,
      simData()$opinion_history[, ncol( simData()$opinion_history)],
      simData()$num_opinions
    )
  })
  
  output$bipartite0Plot <- renderPlot({     
    req(simData())           
    req(input$show_rig0)     
    visual_bipartite(simData()$B0, 
                     simData()$opinion_history[,1], 
                     simData()$num_opinions) })

  output$bipartitePlot <- renderPlot({
    req(simData())          
    req(input$show_rig)     
    visual_bipartite(simData()$B, 
                     simData()$opinion_history[,ncol( simData()$opinion_history)], 
                     simData()$num_opinions) })

  # Voter's dynamics
  output$fracPlot <- renderPlot({
    req(simData())
    visual_step_time(
      simData()$frac_mat,
      simData()$num_opinions
    )
  })

  output$histo <- renderPlot({ 
    req(simData()); 
    visual_histo(simData()$opinion_history, 
                 simData()$num_opinions) })
  
  output$histogramGroup0 <- 
    renderPlot({ req(simData()); 
      visual_histo_pergroup(simData()$opinion_history[,1], 
                            simData()$num_opinions, simData()$members0) })
  
  output$histogramGroup <- renderPlot({ 
    req(simData()); 
    visual_histo_pergroup(simData()$opinion_history[,ncol( simData()$opinion_history)], 
                          simData()$num_opinions, simData()$members) })

  # Epidemic layer
  # Starting reason
  output$startReasonEpi <- renderText({
    req(simData())
    paste("Epidemics started by stability:", simData()$stable)
  })
  
  output$SIR1 <- renderPlot({
    req(simData())
    req(input$runEpidemic)
    visual_step_time_SIR(x = simData()$SIR_df)})

  output$SIR2 <- renderPlot({
    req(simData())
    req(input$runEpidemic)
    visual_step_time_SIR(x = simData()$SIR_df_opinion_red)})

  output$SIR3 <- renderPlot({
    req(simData())
    req(input$runEpidemic)
    visual_step_time_SIR(x = simData()$SIR_df_opinion_blue)})

  output$SIR4 <- renderPlot({
    req(simData())
    req(input$runEpidemic)
    visual_step_time_SIR(x = simData()$SIR_df_grp)})

  output$SIR5 <- renderPlot({
    req(simData())
    req(input$runEpidemic)
    visual_step_time_SIR(x = simData()$SIR_df_out) })

  # Download plots
  output$downloadAllPlots <- downloadHandler(
    filename = function() {
      paste0("all_plots_", Sys.Date(), ".zip")
    },
    content = function(file) {
      req(simData())
      # Temporary folder
      tmpdir <- tempdir()
      
      # Helper function
      save_plot <- function(filename, plot_expr) {
        filepath <- file.path(tmpdir, filename)
        png(
          filename = filepath,
          width = 1400,
          height = 1000,
          res = 150
        )
        plot_expr
        dev.off()
      }

      #==========================
      # Opinion plots
      #==========================

      save_plot(
        "opinion_fractions.png",
        visual_step_time(
          simData()$frac_mat,
          simData()$num_opinions
        )
      )

      save_plot(
        "overall_histogram.png",
        visual_histo(
          simData()$opinion_history,
          simData()$num_opinions
        )
      )

      save_plot(
        "initial_group_histogram.png",
        visual_histo_pergroup(
          simData()$opinion_history[,1],
          simData()$num_opinions,
          simData()$members0
        )
      )

      save_plot(
        "final_group_histogram.png",
        visual_histo_pergroup(
          simData()$opinion_history[, ncol(simData()$opinion_history)],
          simData()$num_opinions,
          simData()$members
        )
      )

      #==========================
      # Epidemic plots
      #==========================

      if (isTRUE(input$runEpidemic)) {

        save_plot(
          "sir_overall.png",
          visual_step_time_SIR(
            x = simData()$SIR_df
          )
        )

        save_plot(
          "sir_red.png",
          visual_step_time_SIR(
            x = simData()$SIR_df_opinion_red
          )
        )

        save_plot(
          "sir_blue.png",
          visual_step_time_SIR(
            x = simData()$SIR_df_opinion_blue
          )
        )

        save_plot(
          "sir_members.png",
          visual_step_time_SIR(
            x = simData()$SIR_df_grp
          )
        )

        save_plot(
          "sir_isolated.png",
          visual_step_time_SIR(
            x = simData()$SIR_df_out
          )
        )
      }

      #==========================
      # Graph plots
      #==========================

      if (isTRUE(input$show_rig0)) {

        save_plot(
          "initial_rig.png",
          visual_step_multi(
            simData()$RIG0,
            simData()$opinion_history[,1],
            simData()$num_opinions
          )
        )

        save_plot(
          "initial_bipartite.png",
          visual_bipartite(
            simData()$B0,
            simData()$opinion_history[,1],
            simData()$num_opinions
          )
        )
      }

      if (isTRUE(input$show_rig)) {

        save_plot(
          "final_rig.png",
          visual_step_multi(
            simData()$RIG,
            simData()$opinion_history[, ncol(simData()$opinion_history)],
            simData()$num_opinions
          )
        )

        save_plot(
          "final_bipartite.png",
          visual_bipartite(
            simData()$B,
            simData()$opinion_history[, ncol(simData()$opinion_history)],
            simData()$num_opinions
          )
        )
      }

      #==========================
      # Create ZIP file
      #==========================

      files <- list.files(
        tmpdir,
        pattern = "\\.png$",
        full.names = TRUE
      )

      zip::zipr(
        zipfile = file,
        files = files
      )
    }
  )



}


#==========================
# --- Run App ---
#==========================
shinyApp(ui, server)

