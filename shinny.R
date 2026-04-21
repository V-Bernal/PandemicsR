library(shiny)
library(igraph)
library(Matrix)
library(PandemicsR)

#==========================
# --- UI ---
#==========================
ui <- fluidPage(
  
  #==========================
  # Title
  #==========================
  titlePanel(h3("Voter - Schelling Multi-Membership Simulation")),
  
  
  #==========================
  # Side panel
  #==========================
  sidebarLayout(
    
    sidebarPanel(
      
      # Section 1: Parameters
      h4("Network parameters"),
      numericInput("n", "Number of individuals", value = 15, min = 5),
      numericInput("m", "Number of groups", value = 3, min = 2),
      numericInput("timesteps", "Gillespie Iterations", value = 1000),
      sliderInput("lambda", "RIG weight parameter lambda", min = 0, step = 0.01, max = 100, value = 10),
      
      # Section 2: Schelling
      h4("Schelling's parameters"),
      checkboxInput("runSchelling", "Schelling model", value = TRUE),
      
      sliderInput("c_param", "Schelling: Edge addition rate param c", min = 0, step = 0.01, max = 1, value = 0.4),
      sliderInput("beta_plus", "Schelling: beta_plus", min = 0, step = 0.01, max = 1, value = 0.5),
      sliderInput("beta_minus", "Schelling: beta_minus", min = 0, step = 0.01, max = 1, value = 0.2),
      sliderInput("T_threshold", "Schelling: T_threshold", min = 0, step = 0.01, max = 1, value = 0.3),
      
      # Section 3: Voter
      h4("Voter's parameters"),
      checkboxInput("runVoter", "Voter model", value = TRUE),
      
      #sliderInput("kappa", "Voter: Poisson rate for opinion update kappa", min = 0, step = 0.01, max = 1, value = 0.3),
      sliderInput("gamma", "Voter: gamma", min = 0, step = 0.01, max = 100, value = 5),
      sliderInput("Numopinions", "Number of Opinions", min = 2, step = 2, max = 4, value = 2),
      
      # Section 4: Extremes
      h4("Extremes"),
      checkboxInput("stubborn_extremes", "Stubbornnes", value = FALSE),
      checkboxInput("deradicalization", "Deradicalization", value = FALSE),
      
      # Section 5: Visualization
      h4("Visualization"),
      checkboxInput("show_rig0", "Show initial Graph", value = FALSE),
      checkboxInput("show_rig", "Show final Graph", value = FALSE),
      
      # Run
      actionButton("runSim", "Run Simulation")
      
          ),
    
    #==========================
    # Main panel: Visualizations
    #==========================
    mainPanel(
      
      tabsetPanel(
        tabPanel("Graphs",
                 plotOutput("rig0Plot", height = "500px"),
                 plotOutput("rigPlot", height = "500px"),
                 plotOutput("bipartite0Plot", height = "500px"),
                 plotOutput("bipartitePlot", height = "500px"),
        ),
        
        tabPanel("Voter dynamics",
                 
                 h3("Overall opinions"),
                 plotOutput("histo"),
                 
                 h3("Time evolution of opinions"),
                 plotOutput("fracPlot", width = "500px", height = "400px")#,
                 #plotOutput("heatmapPlot", height = "350px")
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

#==========================
# --- Server ---
#==========================
server <- function(input, output, session) {

  simData <- eventReactive(input$runSim, {

    #==========================
    # 1. User parameters
    #==========================
    n <- input$n; m <- input$m; t_max <- input$timesteps
    lambda <- input$lambda; c_param <- input$c_param
    gamma <- input$gamma; beta_plus <- input$beta_plus; beta_minus <- input$beta_minus
    T_threshold <- input$T_threshold
    num_opinions <- input$Numopinions

    #==========================
    # 2. Model parameters
    #==========================
    # Homogeneous weights
    ind_w <- rep(lambda, n) #ind_w <- runif(n, 1, 2)
    grp_w <- rep(lambda * n / m, m) #grp_w <- runif(m, 1, 2) * n / m

    # bipartite and RIG
    B0 <- generate_bipartite(n, m, ind_w, grp_w) #B0 <- generate_bipartite(n, m, ind_w, grp_w, lambda)
    RIG0 <- bipartite_to_rig(B0)
    RIG <- RIG0; bipartite <- B0

    # Voter tracker
    opinions <- initialize_opinions_multi(n, num_opinions)
    #opinions0 <- opinions
    opinion_history <- matrix(opinions, ncol = 1)

    levels <- sort(unique(as.vector(opinions)))
    frac_mat <- sapply(levels, function(op) {
      sum(opinions == op)/n
    })
    
    # Group tracker
    members <- vector("list", m)
    red_members <- vector("list", m)
    blue_members <- vector("list", m)
    outsiders <- vector("list", m)
    groups_of_individual <- vector("list", n) # IDs of all groups that individual i belongs to

    for (g in seq_len(m)) {
      ids <- which(bipartite[, g] == 1)
      members[[g]] <- ids
      red_members[[g]] <- ids[opinions[ids] == -1] # moderates
      blue_members[[g]] <- ids[opinions[ids] == +1]
      outsiders[[g]] <- setdiff(seq_len(n), ids)
      for (i in ids) groups_of_individual[[i]] <- c(groups_of_individual[[i]], g)
    }
    # moderates
    Ri <- sapply(red_members, length)
    Bi <- sapply(blue_members, length)
    rig_dirty <- FALSE

    #----- NEW
    event_counter <- 0
    record_every <- n
    members0 <- members
    
    #==========================
    # Combined Voter and Schelling model dynamics Gillespie algorithm
    #==========================
    t <- 0
    while (t < t_max) {

      # The rate of interaction
      Tot <- Ri + Bi
      #voter_term <- gamma * (Ri * Bi/ Tot) # Contact rate ind vd ind: gamma*k
      #voter_term[Tot < 2] <- 0
      
      # new
      voter_term <- numeric(m)
      valid <- Tot >= 2
      voter_term[valid] <- gamma * (Ri[valid] * Bi[valid] / Tot[valid])

      # The rate of joining a group for someone not yet part of a group is c/m.
      # The rate of leaving a group is
      join_term <- (c_param / m) * (n - Tot)/n

      # The rate of leave a group is B+ or B- depending on the thershold
      frac_red <- ifelse(Tot > 0, Ri / Tot, 0)
      frac_blue <- ifelse(Tot > 0, Bi / Tot, 0)
      
      leaveR_rate <- ifelse(frac_red < T_threshold, beta_plus * Ri, beta_minus * Ri)
      leaveB_rate <- ifelse(frac_blue < T_threshold, beta_plus * Bi, beta_minus * Bi)

      # The rate at which the process leaves the state
      lambda_i <- voter_term + join_term + leaveR_rate + leaveB_rate
      lambda_tot <- sum(lambda_i)
      if (lambda_tot <= 0) break

      #====================
      # Group selection at random
      #====================
      group_i <- if(m > 1) sample(1:m, 1L, prob = lambda_i / lambda_tot) else 1

      Ri_g <- Ri[group_i]; Bi_g <- Bi[group_i]
      voter_g <- voter_term[group_i]; join_g <- join_term[group_i]
      leaveR_g <- leaveR_rate[group_i]; leaveB_g <- leaveB_rate[group_i]

      rates_vec <- c(ifelse(voter_g>0, voter_g/2, 0), ifelse(voter_g>0, voter_g/2, 0), join_g, leaveR_g, leaveB_g)
      if (sum(rates_vec)<=0) next

      #================
      # Apply a random move in Gillespie time
      # Gillespie time
      #================
      dt <- rexp(1, lambda_tot)
      t <- t + dt
      if (t >= t_max) break

      #================
      # Move step
      #================
      # 1: Red to Blue. (voter) 
      # 2: Blue to Red. (voter)
      # 3: Add external to group
      # 4: Remove internal Red 
      # 5: Remove internal Blue
      
      #============
      # Trim rate vec according to the number of moves in simulation
      # Voter: moves 1 and 2 
      # Schelling: moves 3 to 5 
      trimmer <- c(rep(input$runVoter, 2), rep(input$runSchelling, 3))
      if (length(trimmer) == length(rates_vec)) { rates_vec <- rates_vec * trimmer } #next
      #============
      
      move <- sample(1:length(rates_vec), 1L, prob = rates_vec / sum(rates_vec))

      if (move==1 && Ri_g>0) { # Red to Blue
        chosen <- sample(red_members[[group_i]],1)
        opinions[chosen] <- +1
        # update every group where individual belongs
        for (g in groups_of_individual[[chosen]]) {
          idx <- which(red_members[[g]]==chosen)
          if (length(idx)>0) {
            red_members[[g]][idx] <- red_members[[g]][length(red_members[[g]])] # Replaces the chosen individual with the last element
            red_members[[g]] <- red_members[[g]][-length(red_members[[g]])] # delete last
            blue_members[[g]] <- c(blue_members[[g]], chosen)
            Ri[g] <- Ri[g]-1
            Bi[g] <- Bi[g]+1
          }
        }

      } else if (move==2 && Bi_g>0) { # Blue to Red
        chosen <- sample(blue_members[[group_i]],1)
        opinions[chosen] <- -1
        # update every group where individual belongs
        for (g in groups_of_individual[[chosen]]) {
          idx <- which(blue_members[[g]]==chosen)
          if (length(idx)>0) {
            blue_members[[g]][idx] <- blue_members[[g]][length(blue_members[[g]])]
            blue_members[[g]] <- blue_members[[g]][-length(blue_members[[g]])]
            red_members[[g]] <- c(red_members[[g]], chosen)
            Bi[g] <- Bi[g]-1
            Ri[g] <- Ri[g]+1
          }
        }

      } else if (move==3 && length(outsiders[[group_i]])>0) { # Join
        chosen_idx <- sample.int(length(outsiders[[group_i]]),1)
        chosen <- outsiders[[group_i]][chosen_idx]
        outsiders[[group_i]][chosen_idx] <- outsiders[[group_i]][length(outsiders[[group_i]])]
        outsiders[[group_i]] <- outsiders[[group_i]][-length(outsiders[[group_i]])]

        members[[group_i]] <- c(members[[group_i]], chosen)
        groups_of_individual[[chosen]] <- c(groups_of_individual[[chosen]], group_i)
        if (opinions[chosen]==-1) {
          red_members[[group_i]] <- c(red_members[[group_i]], chosen)
          Ri[group_i] <- Ri[group_i]+1
        } else {
          blue_members[[group_i]] <- c(blue_members[[group_i]], chosen)
          Bi[group_i] <- Bi[group_i]+1
        }
        rig_dirty <- TRUE

      } else if (move==4 && Ri_g>0) {
        chosen_idx <- sample.int(Ri_g,1)
        chosen <- red_members[[group_i]][chosen_idx]
        red_members[[group_i]][chosen_idx] <- red_members[[group_i]][length(red_members[[group_i]])]
        red_members[[group_i]] <- red_members[[group_i]][-length(red_members[[group_i]])]

        idx_m <- which(members[[group_i]]==chosen)
        members[[group_i]][idx_m] <- members[[group_i]][length(members[[group_i]])]
        members[[group_i]] <- members[[group_i]][-length(members[[group_i]])]
        groups_of_individual[[chosen]] <- setdiff(groups_of_individual[[chosen]], group_i)
        outsiders[[group_i]] <- c(outsiders[[group_i]], chosen)
        Ri[group_i] <- Ri[group_i]-1
        rig_dirty <- TRUE

      } else if (move==5 && Bi_g>0) {
        chosen_idx <- sample.int(Bi_g,1)
        chosen <- blue_members[[group_i]][chosen_idx]
        blue_members[[group_i]][chosen_idx] <- blue_members[[group_i]][length(blue_members[[group_i]])]
        blue_members[[group_i]] <- blue_members[[group_i]][-length(blue_members[[group_i]])]

        idx_m <- which(members[[group_i]]==chosen)
        members[[group_i]][idx_m] <- members[[group_i]][length(members[[group_i]])]
        members[[group_i]] <- members[[group_i]][-length(members[[group_i]])]
        groups_of_individual[[chosen]] <- setdiff(groups_of_individual[[chosen]], group_i)
        outsiders[[group_i]] <- c(outsiders[[group_i]], chosen)
        Bi[group_i] <- Bi[group_i]-1
        rig_dirty <- TRUE
      }

      # ---- record step ----
      event_counter <- event_counter + 1

      #if (event_counter %% record_every == 0) {
      #  opinion_history <- cbind(opinion_history, opinions)
      #}
      
      if (event_counter %% record_every == 0) {
        frac_temp <- sapply(levels, function(op) {
          sum(opinions == op)/n
        })
        frac_mat <- cbind(frac_mat, frac_temp)
      }
      
    }

    # Final opinion history
    opinion_history <- cbind(opinion_history, opinions)
    frac_mat <- t(frac_mat)
    colnames(frac_mat) <- levels
    
    # if (rig_dirty) {
    #   bipartite <- reconstruct_bipartite(members, n, m)
    #   RIG <- bipartite_to_rig(bipartite)
    # }

    list(
      
      # comment
      B0=B0, B=bipartite, RIG=RIG, RIG0=RIG0,
      opinions=opinions, #opinions0=opinions0,
      #
      members0=members0,frac_mat=frac_mat,
      opinion_history=opinion_history,
      num_opinions=num_opinions,members=members
    )
  })

  #==========================
  # --- Output ---
  #==========================
  
  #---------------
  # comment
  #output$rig0Plot <- renderPlot({ req(simData()) ; visual_step_multi(simData()$RIG0, simData()$opinion_history[,1], simData()$num_opinions) })
  
  output$rig0Plot <- renderPlot({
    req(simData())           # simulation must exist
    req(input$show_rig0)     # checkbox must be TRUE
    
    visual_step_multi(
      simData()$RIG0,
      simData()$opinion_history[,1],
      simData()$num_opinions
    )
  })
  
  output$rigPlot <- renderPlot({
    req(simData())           # simulation must exist
    req(input$show_rig)     # checkbox must be TRUE
    
    visual_step_multi(
      simData()$RIG,
      simData()$opinion_history[, ncol( simData()$opinion_history)],
      simData()$num_opinions
    )
  })
  
  # output$rigPlot <- renderPlot({ req(simData()); visual_step_multi(simData()$RIG,
  #                                                                  simData()$opinion_history[,ncol( simData()$opinion_history)], simData()$num_opinions) })
  output$bipartite0Plot <- renderPlot({     req(simData())           # simulation must exist
    req(input$show_rig0)     # checkbox must be TRUE; 
    visual_bipartite(simData()$B0, simData()$opinion_history[,1], simData()$num_opinions) })
  
  output$bipartitePlot <- renderPlot({     
    req(simData())           # simulation must exist
    req(input$show_rig)     # checkbox must be TRUE; 
    visual_bipartite(simData()$B, simData()$opinion_history[,ncol( simData()$opinion_history)], simData()$num_opinions) })
  
  #---------------
  output$fracPlot <- renderPlot({ req(simData()); visual_step_time(simData()$frac_mat, simData()$num_opinions, length(simData()$opinion_history[,1])) })
  output$histo <- renderPlot({ req(simData()); visual_histo(simData()$opinion_history, simData()$num_opinions) })
  #output$heatmapPlot <- renderPlot({ req(simData()); heatmapPlot(simData()$opinion_history, simData()$num_opinions) })
  #output$histogramGroup0 <- renderPlot({ req(simData()); visual_histo_pergroup(simData()$opinion_history[,1], simData()$num_opinions, simData()$B0) })
  #output$histogramGroup <- renderPlot({ req(simData()); visual_histo_pergroup(simData()$opinion_history[,ncol( simData()$opinion_history)], simData()$num_opinions, simData()$B) })
  output$histogramGroup0 <- renderPlot({ req(simData()); visual_histo_pergroup(simData()$opinion_history[,1], simData()$num_opinions, simData()$members0) })
  output$histogramGroup <- renderPlot({ req(simData()); visual_histo_pergroup(simData()$opinion_history[,ncol( simData()$opinion_history)], simData()$num_opinions, simData()$members) })
  
}

#==========================
# --- Run App ---
#==========================
shinyApp(ui, server)
