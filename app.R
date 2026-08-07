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
        
        fluidRow(column(6, numericInput("timesteps", "Max Gillespie time", value = 100)),
          column(6,numericInput("lambda", "RIG weight (lambda) ",value = 1)) )
              ),

      # Section 2: Schelling
      wellPanel(
        h4("Schelling model"),
        checkboxInput("runSchelling", "Schelling model", value = FALSE),
        
        conditionalPanel(
          condition = "input.runSchelling",
          
        sliderInput("c_param", "join group rate c", min = 0, step = 0.01, max = 1, value = 0.4),

        h6("*to scale the join rate: tick the box"),
        checkboxInput("scaled_n", "scaled by N", value = FALSE),
        checkboxInput("scaled_m", "scaled by m", value = FALSE),

        sliderInput("beta_plus", "leave rate above tolerance", min = 0, step = 0.01, max = 1, value = 0.5),
        sliderInput("beta_minus", "leave rate below tolerance", min = 0, step = 0.01, max = 1, value = 0.2),
        sliderInput("T_threshold", "tolerance threshold", min = 0, step = 0.01, max = 1, value = 0.3)        )
        ),

      # Section 3: Voter
      wellPanel(
        h4("Voter's model"),
        checkboxInput("runVoter", "Voter model", value = FALSE),
        
        conditionalPanel(
          condition = "input.runVoter",
        #sliderInput("kappa", "Voter: Poisson rate for opinion update kappa", min = 0, step = 0.01, max = 1, value = 0.3),
        sliderInput("gamma", "Opinions rate", min = 0, step = 0.1, max = 100, value = 5),
        sliderInput("Numopinions", "Number of Opinions", min = 2, step = 2, max = 4, value = 4),

        # Section 4: Extremes
        #h5("Extremes"),
        sliderInput("alpha", "radicalization rate", min = 0, step = 0.01, max = 1, value = 0.1),
        sliderInput("alpha_deradicalization", "deradicalization rate*", min = 0, step = 0.01, max = 1, value = 0),
        h6("* for stubborn opinions: set de-radicalization to zero"),
        sliderInput("alpha0", "spontaneous (de)-radicalization rate", min = 0, step = 0.01, max = 1, value = 0)
        )
        ),

      # Section 5: Epidemics
      wellPanel(
        h4("Epidemics model"),
        checkboxInput("runEpidemic", "Epidemic model", value = FALSE),
        
        conditionalPanel(
          
        condition = "input.runEpidemic",
        sliderInput("I0", "Fraction of Infected", min = 0.01, step = 0.1, max = 1, value = 0.1),
        sliderInput("gamma_epi", "Recovery rate", min = 0, step = 0.01, max = 1, value = 0.3),
        sliderInput("beta_red", "Infection rate red", min = 0, step = 0.01, max = 1, value = 0.5),
        sliderInput("beta_blue", "Infection rate blue", min = 0, step = 0.01, max = 1, value = 0.2)
        
        ),
      
      # Section 5.1: Epidemics activation
      wellPanel(
        h4("Epidemics activation"),
      
      selectInput(
        "epi_trigger",
        "Start epidemic after:",
        choices = c(
          "Fixed time" = "time",
          "Fixed events" = "events",
          "Events per node" = "events_per_node",
          "Opinion stability" = "opinion_stability",
          "Segregation stability" = "segregation_stability"
        ),
        selected = "time"
      ),
      
      conditionalPanel(
        condition = "input.epi_trigger == 'time'",
        numericInput(
          "epi_time",
          "Epidemic start time:",
          value = 0,
          min = 0
        )
      ),
      
      conditionalPanel(
        condition = "input.epi_trigger == 'events'",
        numericInput(
          "epi_events",
          "Number of warm-up events:",
          value = 0,
          min = 0, max = 0
        )
      ),
      
      conditionalPanel(
        condition = "input.epi_trigger == 'events_per_node'",
        numericInput(
          "epi_events_node",
          "Warm-up events per node:",
          value = 0,
          min = 1, max = 1
        )
      ),
      
      conditionalPanel(
        condition = "input.epi_trigger == 'opinion_stability'",
        numericInput(
          "epi_epsilon_opinion",
          "Opinion change threshold:",
          value = 0,
          min = 0, max = 0
        )
      ),
      
      conditionalPanel(
        condition = "input.epi_trigger == 'segregation_stability'",
        numericInput(
          "epi_epsilon_seg",
          "Segregation change threshold:",
          value = 0,
          min = 0, max = 0
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


    # params <- list(
    #   n = 15,
    #   m = 3,
    #   t_max = 1000,
    #   lambda = 1,
    #   gamma = 5,
    #   beta_plus = 0.2,
    #   beta_minus = 0.2,
    #   T_threshold = 1000,
    #   num_opinions = 2,
    #   c_param = 0.4,
    #   alpha = 0.5,
    #   alpha_deradicalization = 0.5,
    #   alpha0 = 0.5,
    #   runVoter = TRUE,
    #   runSchelling = TRUE,
    #   runEpidemic = FALSE,
    #   beta_red = 0.2,
    #   beta_blue = 0.2,
    #   gamma_epi = 0.5,
    #   I0 = 0.5
    # )

    params <- list(
      n = input$n,
      m = input$m,
      t_max = input$timesteps,
      lambda = input$lambda,
      gamma = input$gamma,
      beta_plus = input$beta_plus,
      beta_minus = input$beta_minus,
      c_param = input$c_param,
      T_threshold = input$T_threshold,
      num_opinions = input$Numopinions,
      alpha = input$alpha,
      alpha_deradicalization = input$alpha_deradicalization,
      alpha0 = input$alpha0,
      runEpidemic = input$runEpidemic,
      runVoter = input$runVoter,
      runSchelling = input$runSchelling,
      beta_red = input$beta_red,
      beta_blue = input$beta_blue,
      gamma_epi = input$gamma_epi,
      I0 = input$I0
    )
    
    #params$epi_trigger <- input$epi_trigger[input$epi_trigger == c("time", "events", "events_per_node", "opinion_stability", "segregation_stability")]
    params$epi_time <- input$epi_time
    #params$epi_events <- input$epi_events
    
    # Simulation
    run_simulation(params)

    }

    )


  #==========================
  # Section 6: Visualization
  #==========================
  
  # Stopping reason
  output$stopReason <- renderText({
    req(simData())
    paste("Simulation stopped:", simData()$stop_reason,
          "| Final Gillespie time:", round(simData()$final_time, 2))
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

