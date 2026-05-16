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
        numericInput("n", "Number of individuals", value = 15, min = 5),
        numericInput("m", "Number of groups", value = 3, min = 2),
        numericInput("timesteps", "Gillespie Iterations", value = 1000),
        numericInput("lambda", "RIG weight parameter lambda", value = 1, min = 0)
      ),

      # Section 2: Schelling
      wellPanel(
        h4("Schelling model"),
        checkboxInput("runSchelling", "Schelling model", value = FALSE),

        sliderInput("c_param", "Schelling: Edge addition rate param c", min = 0, step = 0.01, max = 1, value = 0.4),
        sliderInput("beta_plus", "Schelling: beta_plus", min = 0, step = 0.01, max = 1, value = 0.5),
        sliderInput("beta_minus", "Schelling: beta_minus", min = 0, step = 0.01, max = 1, value = 0.2),
        sliderInput("T_threshold", "Schelling: T_threshold", min = 0, step = 0.01, max = 1, value = 0.3)
      ),

      # Section 3: Voter
      wellPanel(
        h4("Voter's model"),
        checkboxInput("runVoter", "Voter model", value = FALSE),
        #sliderInput("kappa", "Voter: Poisson rate for opinion update kappa", min = 0, step = 0.01, max = 1, value = 0.3),
        sliderInput("gamma", "Voter: gamma", min = 0, step = 0.01, max = 100, value = 5),
        sliderInput("Numopinions", "Number of Opinions", min = 2, step = 2, max = 4, value = 4),

        # Section 4: Extremes
        #h5("Extremes"),
        sliderInput("alpha", "radicalization rate", min = 0, step = 0.01, max = 1, value = 0),
        sliderInput("alpha_deradicalization", "deradicalization rate*", min = 0, step = 0.01, max = 1, value = 0),
        h6("*Stubbornnes: set deradicalization to zero")
      ),

      # Section 5: Epidemics
      wellPanel(
        h4("Epidemics model"),
        checkboxInput("runEpidemic", "Epidemic model", value = FALSE),
        sliderInput("I0", "Number of Infected", min = 1, step = 1, max = 20, value = 5),
        sliderInput("gamma_epi", "Epidemic: recovery rate", min = 0, step = 0.01, max = 1, value = 0.3),
        sliderInput("beta_red", "Epidemic: infection rate red", min = 0, step = 0.01, max = 1, value = 0.5),
        sliderInput("beta_blue", "Epidemic: infection rate blue", min = 0, step = 0.01, max = 1, value = 0.2)
      ),

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

    # real time
    start_time <- Sys.time()

    #==========================
    # 1. User parameters
    #==========================
    n <- input$n; m <- input$m; t_max <- input$timesteps
    lambda <- input$lambda; c_param <- input$c_param
    gamma <- input$gamma; beta_plus <- input$beta_plus; beta_minus <- input$beta_minus
    T_threshold <- input$T_threshold
    num_opinions <- input$Numopinions

    # radicalization
    alpha <- input$alpha
    alpha_deradicalization <- input$alpha_deradicalization

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
    
    #==========================
    # Model trackers
    #==========================
    # Voter tracker
    opinions <- initialize_opinions_multi(n, num_opinions)
    opinion_history <- matrix(opinions, ncol = 1)

    #is_member <- logical(n)
    #membership_count <- integer(n)

    levels <- if(num_opinions == 4) c(-2, -1, 1, 2) else c(-1, 1)# levels <- sort(unique(as.vector(opinions)))
    
    frac_mat <- sapply(levels, function(op) {
      sum(opinions == op)/n
    })

    # Group tracker moderates
    members <- vector("list", m)
    red_members <- vector("list", m)
    blue_members <- vector("list", m)
    outsiders <- vector("list", m)
    groups_of_individual <- vector("list", n) # IDs of all groups that individual i belongs to
    
    # Moderates
    for (g in seq_len(m)) {
      #group
      ids <- which(bipartite[, g] == 1)
      members[[g]] <- ids
      #opinion in group
      red_members[[g]] <- ids[opinions[ids] == -1]
      blue_members[[g]] <- ids[opinions[ids] == +1]
      # outside the group is not necessarily isolated
      outsiders[[g]] <- setdiff(seq_len(n), ids)
      for (i in ids) groups_of_individual[[i]] <- c(groups_of_individual[[i]], g)
    }
    
    # moderates
    Ri <- sapply(red_members, length)
    Bi <- sapply(blue_members, length)
    rig_dirty <- FALSE

    #----- NEW
    event_counter <- 0
    record_every <- 1#floor(1e4 / 8 / (16 + num_opinions))#max(1, floor( t_max/ 10))
    members0 <- members

    # group members
    all_group_members <- unique(unlist(members))
    #all_group_members <- which(membership_count > 0)
    
    # Isolated Outsiders all
    outsiders_all <- setdiff(seq_len(n), all_group_members)

    
    #--time of infection----------
    # Initial time
    time_hist <- c(0)
    
    inf_time <- c()
    inf_camp <- c()
    
    stop_reason <- "Reached maximum time"
    
    #=====================
    # Initialize epidemic (O(1) structure)
    #=====================
    if (isTRUE(input$runEpidemic)) {
    #==========================
    # New Epidemic states
    #==========================
    # Each ind can be in one state S or I or R
    S <- 0; I <- 1; R <- 2
    
      total_red  <- sum(opinions < 0)
      total_blue <- sum(opinions > 0)
    
      # Initial Parameters (NEW)
      I0 <- input$I0 # number of initial infected
      beta_red  <- input$beta_red
      beta_blue <- input$beta_blue # infection rate blue opinion
      gamma_epi <- input$gamma_epi #recovery rate
    
      epi <- rep(S, n)
      initial_infected <- sample.int(n, I0)
      epi[initial_infected] <- I
      
      # Maintain explicit sets
      S_nodes <- setdiff(seq_len(n), initial_infected)
      I_nodes <- initial_infected
      
      S_count <- length(S_nodes)
      I_count <- length(I_nodes)
      R_count <- 0
      
      # Positions for O(1) removal
      pos_in_S <- integer(n)
      pos_in_I <- integer(n)
      
      for (k in seq_along(S_nodes)) pos_in_S[S_nodes[k]] <- k
      for (k in seq_along(I_nodes)) pos_in_I[I_nodes[k]] <- k
      
      #=====================
      # Initialize epidemic history trackers
      #=====================
      
      # Overall SIR counts
      S_hist <- c(S_count)
      I_hist <- c(I_count)
      R_hist <- c(R_count)
      
      #=====================
      # Opinion camp SIR
      #=====================
      
      # Initial recovered counts
      R_red  <- 0
      R_blue <- 0
      
      # Initial infected counts
      I_red  <- sum(opinions[I_nodes] < 0)
      I_blue <- sum(opinions[I_nodes] > 0)
      
      # Initial susceptible counts
      S_red  <- sum(opinions < 0) - I_red
      S_blue <- sum(opinions > 0) - I_blue
      
      # Store history
      S_hist_red  <- c(S_red)
      I_hist_red  <- c(I_red)
      R_hist_red  <- c(R_red)
      
      S_hist_blue <- c(S_blue)
      I_hist_blue <- c(I_blue)
      R_hist_blue <- c(R_blue)
      
      #=====================
      # Group membership SIR
      #=====================

      
      S_hist_grp <- c(sum(epi[all_group_members] == S))
      I_hist_grp <- c(sum(epi[all_group_members] == I))
      R_hist_grp <- c(sum(epi[all_group_members] == R))
      
      S_hist_out <- c(sum(epi[outsiders_all] == S))
      I_hist_out <- c(sum(epi[outsiders_all] == I))
      R_hist_out <- c(sum(epi[outsiders_all] == R))
      
    }


    #==========================
    # Gillespie algorithm
    #==========================
    t <- 0
    while (t < t_max) {

      # stopping time
      if (isTRUE(input$runEpidemic) && all(epi == R)) {
        stop_reason <- "All individuals recovered (epidemic ended)"
        break
      }

      if (num_opinions == 4 && all(opinions %in% c(-2, 2))) {
        stop_reason <- "Full polarization (all extreme opinions)"
        break
      }
      
      #==========================
      # Voter's rate of interaction among moderates
      #==========================
      Tot <- Ri + Bi
      voter_term <- numeric(m)
      valid <- (Ri > 0 & Bi > 0) #Tot >= 2
      voter_term[valid] <- gamma * (Ri[valid] * Bi[valid] / Tot[valid])

      #==========================
      # joining a group rates
      ##==========================
      # joining a group for someone not yet part of a group is c/m.
      join_term <- (c_param / m) * sapply(outsiders, length)/n #(n - Tot)/n

      #==========================
      # leaving a group rates
      ##==========================
      # The rate of leaving a group is B+ or B- depending on the threshold
      frac_red <- ifelse(Tot > 0, Ri / Tot, 0)
      frac_blue <- ifelse(Tot > 0, Bi / Tot, 0)

      leaveR_rate <- ifelse(frac_red < T_threshold, beta_plus * Ri, beta_minus * Ri)
      leaveB_rate <- ifelse(frac_blue < T_threshold, beta_plus * Bi, beta_minus * Bi)

      #==========================
      # Epidemic rates (NEW)
      #==========================
      infection_rate_tot <- 0
      recovery_rate_tot <- 0
      
      if (isTRUE(input$runEpidemic)) {
        
        #I_count <- sum(epi == I)
        prevalence <- I_count / n # global = random mixing
        
        # Opinion camp-based infection rates
        beta_vec <- ifelse(opinions < 0, beta_red, beta_blue)
        
        # Infection: only for S individuals
        ##infection_rates <- ifelse(epi == S, beta_vec * prevalence, 0)
        
        # Recovery: only for I individuals
        ##recovery_rates <- ifelse(epi == I, gamma_epi, 0)
        
        # Aggregate epidemic rates
        #infection_rate_tot <- sum(infection_rates)
        #recovery_rate_tot  <- sum(recovery_rates)
        
        infection_rate_tot <- (I_count / n) * sum(beta_vec[S_nodes])
        recovery_rate_tot  <- gamma_epi * I_count
        
        } 
      
      #==========================
      # New: de-radicalize rates
      #==========================
      rate_radicalize_red  <- numeric(m)
      rate_radicalize_blue <- numeric(m)
      rate_deradicalize_red <- numeric(m)
      rate_deradicalize_blue <- numeric(m)


      for (g in seq_len(m)) {
        num_extreme_red   <- sum(opinions[members[[g]]] == -2)
        num_extreme_blue  <- sum(opinions[members[[g]]] ==  2)
        num_moderate_red  <- sum(opinions[members[[g]]] == -1)
        num_moderate_blue <- sum(opinions[members[[g]]] ==  1)

        Tot_g <- length(members[[g]])
        factor <- if (Tot_g > 0) 1 / Tot_g else 0 # avoid NAN or division by zero

        rate_radicalize_red[g]  <- alpha * num_extreme_red  * num_moderate_red  * factor
        rate_radicalize_blue[g] <- alpha * num_extreme_blue * num_moderate_blue * factor

        rate_deradicalize_red[g]  <- alpha_deradicalization * num_extreme_red  * num_moderate_red  * factor
        rate_deradicalize_blue[g] <- alpha_deradicalization * num_extreme_blue * num_moderate_blue * factor

      }

      #==========================
      # Global state rate
      #==========================
      # Each group has a lambda state rate
      lambda_i <- voter_term + join_term + leaveR_rate + leaveB_rate +
        rate_radicalize_red + rate_radicalize_blue +
        rate_deradicalize_red + rate_deradicalize_blue

      social_rate <- sum(lambda_i)
      
      epi_rate <- infection_rate_tot + recovery_rate_tot

      lambda_tot <- social_rate + epi_rate

      if (lambda_tot <= 0) break
      
      #============
      # Trim rate vec according to the number of moves in simulation
      enabled_moves <- c()
      
      if (isTRUE(input$runVoter)) {
        enabled_moves <- c(enabled_moves, 1, 2)
      }
      
      if (isTRUE(input$runSchelling)) {
        enabled_moves <- c(enabled_moves, 3, 4, 5)
      }
      
      if (!is.null(num_opinions) && num_opinions == 4) {
        enabled_moves <- c(enabled_moves, 6, 7)
      }
      
      if (alpha_deradicalization > 0) {
        enabled_moves <- c(enabled_moves, 8, 9)
      }
      
      #================
      # Apply a random move in Gillespie time
      # Gillespie time
      #================
      dt <- rexp(1, lambda_tot)
      t <- t + dt
      if (t >= t_max) break

      #====================
      # Select Social or Epidemic event
      #====================
      u_event <- runif(1)

      if (u_event < social_rate / lambda_tot) {

        #===========
        # Social event
        #===========
        
        #====================
        # Group g selection at random
        #====================
        group_i <- NULL
        
        if (social_rate > 0) {
          group_i <- if (m > 1) {
            sample(1:m, 1L, prob = lambda_i / social_rate)
          } else 1
        }
        
        # subset the moderate rates per group
        Ri_g <- Ri[group_i]; Bi_g <- Bi[group_i]
        voter_g <- voter_term[group_i]; join_g <- join_term[group_i]
        leaveR_g <- leaveR_rate[group_i]; leaveB_g <- leaveB_rate[group_i]
        
        rates_vec <- c(ifelse(voter_g>0, voter_g/2, 0),
                       ifelse(voter_g>0, voter_g/2, 0),
                       join_g, leaveR_g, leaveB_g)
        
        # subset the extreme rates per group
        rates_vec <- c(rates_vec,
                       rate_radicalize_red[group_i], # move 6
                       rate_radicalize_blue[group_i] # move 7
        )
        
        rates_vec <- c(rates_vec,
                       rate_deradicalize_red[group_i], # move 8
                       rate_deradicalize_blue[group_i] # move 9
        )
        
        if (sum(rates_vec)<=0) next
        
        
        # Filter rates
        rates_sub <- rates_vec[enabled_moves]
        
        # Safety checks
        if (length(rates_sub) == 0 || sum(rates_sub) <= 0) next
        
        #====================
        # Move selection at random
        #====================
        # 1: Red to Blue. (voter)
        # 2: Blue to Red. (voter)
        # 3: Add external to group
        # 4: Remove internal Red
        # 5: Remove internal Blue
        #============
        # 6 and 7: moderate turns extreme of same color (voter)
        # 8 and 9:  extreme to moderate of same color (voter)
        #============
        move <- sample(enabled_moves, 1L, prob = rates_sub / sum(rates_sub))


        if (move == 1 && (Ri_g > 0 && Bi_g > 0)) {
          # update Red to Blue
          chosen_idx <- sample.int(length(red_members[[group_i]]), 1)
          chosen <- red_members[[group_i]][chosen_idx]
          #chosen <- sample(red_members[[group_i]],1)
          opinions[chosen] <- +1

          # update epidemic camp counters
          if (epi[chosen] == I) {
          
             I_red  <- I_red - 1
             I_blue <- I_blue + 1
          
           } else if (epi[chosen] == R) {
          
             R_red  <- R_red - 1
             R_blue <- R_blue + 1
           }

          # global opinion counters
          total_red  <- total_red - 1
          total_blue <- total_blue + 1

          # derived susceptible counters
          S_red  <- total_red  - I_red - R_red
          S_blue <- total_blue - I_blue - R_blue

          # update opinion in every group
          for (g in groups_of_individual[[chosen]]) {
            idx <- which(red_members[[g]] == chosen)
            if (length(idx)>0) {
              red_members[[g]][idx] <- red_members[[g]][length(red_members[[g]])] # Replaces the chosen individual with the last element
              red_members[[g]] <- red_members[[g]][-length(red_members[[g]])] # delete last
              blue_members[[g]] <- c(blue_members[[g]], chosen)
              Ri[g] <- Ri[g]-1
              Bi[g] <- Bi[g]+1
            }
          }



        } else if (move == 2 && (Ri_g > 0 && Bi_g > 0)) {

          # update Blue to Red
          chosen_idx <- sample.int(length(blue_members[[group_i]]), 1)
          chosen <- blue_members[[group_i]][chosen_idx]
          # chosen <- sample(blue_members[[group_i]],1)
          opinions[chosen] <- -1

          # update epidemic camp counters
          if (epi[chosen] == I) {

            I_blue <- I_blue - 1
            I_red  <- I_red + 1

          } else if (epi[chosen] == R) {

            R_blue <- R_blue - 1
            R_red  <- R_red + 1
          }

          # global opinion counters
          total_red  <- total_red + 1
          total_blue <- total_blue - 1

          # derived susceptible counters
          S_red  <- total_red  - I_red - R_red
          S_blue <- total_blue - I_blue - R_blue

          # update in every group
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

        } else if (move == 3 && (length( outsiders[[group_i]] ) > 0) ) {
          
          # Join external to group
          chosen_idx <- sample.int(length(outsiders[[group_i]]),1)
          chosen <- outsiders[[group_i]][chosen_idx]
          outsiders[[group_i]][chosen_idx] <- outsiders[[group_i]][length(outsiders[[group_i]])]
          outsiders[[group_i]] <- outsiders[[group_i]][-length(outsiders[[group_i]])]

          members[[group_i]] <- c(members[[group_i]], chosen)
          groups_of_individual[[chosen]] <- c(groups_of_individual[[chosen]], group_i)


          if (opinions[chosen] == -1) {
            red_members[[group_i]] <- c(red_members[[group_i]], chosen)
            Ri[group_i] <- Ri[group_i] + 1
          } else if (opinions[chosen] == 1){
            blue_members[[group_i]] <- c(blue_members[[group_i]], chosen)
            Bi[group_i] <- Bi[group_i] + 1
          }

          rig_dirty <- TRUE

        } else if (move == 4 && Ri_g > 0) {
          # Remove internal moderate Red
          chosen_idx <- sample.int(length(red_members[[group_i]]), 1)
          chosen <- red_members[[group_i]][chosen_idx]
          red_members[[group_i]][chosen_idx] <- red_members[[group_i]][length(red_members[[group_i]])]
          red_members[[group_i]] <- red_members[[group_i]][-length(red_members[[group_i]])]

          idx_m <- which(members[[group_i]] == chosen)
          members[[group_i]][idx_m] <- members[[group_i]][length(members[[group_i]])]
          members[[group_i]] <- members[[group_i]][-length(members[[group_i]])]
          
          groups_of_individual[[chosen]] <- setdiff(groups_of_individual[[chosen]], group_i)
          outsiders[[group_i]] <- c(outsiders[[group_i]], chosen)
          
          Ri[group_i] <- Ri[group_i] - 1
          rig_dirty <- TRUE

        } else if (move == 5 && Bi_g > 0) {
          # Remove internal moderate blue
          chosen_idx <- sample.int(length(blue_members[[group_i]]),1)
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

          } else if (move == 6 &&
                   rate_radicalize_red[group_i] > 0 &&
                   Ri[group_i] > 0) {
            # moderate to extreme red
          moderates <- members[[group_i]][opinions[members[[group_i]]] == -1]
          extremes <- members[[group_i]][opinions[members[[group_i]]] == -2]

          if (length(moderates) > 0 && length(extremes) > 0) { #check if we need length(extremes) > 0
            
            chosen_moderate <- if(length(moderates) == 1) moderates else sample(moderates, 1)
            #chosen_moderate <- sample(moderates, 1)
            opinions[chosen_moderate] <- -2

            # if (epi[chosen_moderate] == I) {
            #   # stays red camp → NO counter change
            # }
          
          
          # UPDATE ALL GROUPS OF THE INDIVIDUAL
          for (g in groups_of_individual[[chosen_moderate]]) {
            idx <- which(red_members[[g]] == chosen_moderate)
            if (length(idx) > 0) {
              red_members[[g]][idx] <- red_members[[g]][length(red_members[[g]])]
              red_members[[g]] <- red_members[[g]][-length(red_members[[g]])]
              Ri[g] <- Ri[g] - 1 # moderate red turns extreme red
            }
          }
          
          }
          
          } else if (move == 7 &&
                   rate_radicalize_blue[group_i] > 0 &&
                   Bi[group_i] > 0) {
          # moderate to extreme blue

          moderates <- members[[group_i]][opinions[members[[group_i]]] == 1]
          extremes <- members[[group_i]][opinions[members[[group_i]]] == 2]
          
          if (length(moderates) > 0 && length(extremes) > 0) { #check if we need length(extremes) > 0
            
            chosen_moderate <- if(length(moderates) == 1) moderates else sample(moderates, 1)
            #chosen_moderate <- sample(moderates, 1)
            opinions[chosen_moderate] <- 2
          #}
          
          # UPDATE ALL GROUPS OF THE INDIVIDUAL
          for (g in groups_of_individual[[chosen_moderate]]) {
            idx <- which(blue_members[[g]] == chosen_moderate)
            if (length(idx) > 0) {
              blue_members[[g]][idx] <- blue_members[[g]][length(blue_members[[g]])]
              blue_members[[g]] <- blue_members[[g]][-length(blue_members[[g]])]
              Bi[g] <- Bi[g] - 1 # moderate blue turns extreme blues
            }
          }
          }  
          
        }  else if (move == 8 &&
                   rate_deradicalize_red[group_i] > 0 &&
                   Ri[group_i] > 0) {
          # extreme to moderate red

          moderates <- members[[group_i]][opinions[members[[group_i]]] == -1]
          extremes <- members[[group_i]][opinions[members[[group_i]]] == -2]

          if (length(moderates) > 0 && length(extremes) > 0) {

            #chosen_extreme <- sample(extremes, 1)
            chosen_extreme <- if(length(extremes) == 1) extremes else sample(extremes, 1)
            opinions[chosen_extreme] <- -1

            # if (epi[chosen_extreme] == I) {
            #   # stays  → NO counter change
            # }
            
            # UPDATE ALL GROUPS OF THE INDIVIDUAL
            for (g in groups_of_individual[[chosen_extreme]]) {
              
              # add to red_members (moderates)
              red_members[[g]] <- c(red_members[[g]], chosen_extreme)
              
              # update count
              Ri[g] <- Ri[g] + 1
            }
          }

          }  else if (move == 9 &&
                   rate_deradicalize_blue[group_i] > 0 &&
                   Bi[group_i] > 0) {
          # extreme to moderate blue

          moderates <- members[[group_i]][opinions[members[[group_i]]] == 1]
          extremes <- members[[group_i]][opinions[members[[group_i]]] == 2]

          if (length(moderates) > 0 && length(extremes) > 0) {

            chosen_extreme <- if(length(extremes) == 1) extremes else sample(extremes, 1)
            #chosen_extreme <- sample(extremes, 1)
            opinions[chosen_extreme] <- 1

            # if (epi[chosen_extreme] == I) {
            #   # stays  → NO counter change
            # }
            
            # UPDATE ALL GROUPS OF THE INDIVIDUAL
            for (g in groups_of_individual[[chosen_extreme]]) {
              
              # add to blue_members (moderates)
              blue_members[[g]] <- c(blue_members[[g]], chosen_extreme)
              
              # update count
              Bi[g] <- Bi[g] + 1
            }
            
          }
        
      }
      
      } else {

        #==========================
        # Section 5. Epidemic event
        #==========================

        if ((infection_rate_tot + recovery_rate_tot) == 0) next

        # Decide infection vs recovery
        if (runif(1) < infection_rate_tot / (infection_rate_tot + recovery_rate_tot)) {

          # # # -------- Infection event --------
          # pick random susceptible
          if (length(S_nodes) == 0) next
          weights <- beta_vec[S_nodes]

          if (sum(weights) <= 0) next
          k <- sample.int(length(S_nodes), 1, prob = weights)

          i <- S_nodes[k]
          stopifnot(epi[i] == S)

          # opinion at time of infection
          camp_i <- ifelse(opinions[i] < 0, "red", "blue")
          inf_time <- c(inf_time, t)
          inf_camp <- c(inf_camp, camp_i)

          # remove from S (swap-delete)
          # remove from S_nodes
          last <- S_nodes[length(S_nodes)]

          S_nodes[k] <- last
          pos_in_S[last] <- k

          S_nodes <- S_nodes[-length(S_nodes)]
          pos_in_S[i] <- 0

          # add to I_nodes
          I_nodes <- c(I_nodes, i)
          pos_in_I[i] <- length(I_nodes)

          epi[i] <- I

          S_count <- S_count - 1
          I_count <- I_count + 1

          if (opinions[i] < 0) {
            I_red <- I_red + 1
          } else {
            I_blue <- I_blue + 1
          }

          S_red  <- total_red  - I_red - R_red
          S_blue <- total_blue - I_blue - R_blue

        } else {

          # -------- Recovery event --------
          # pick random infected
          if (length(I_nodes) == 0) next
          k <- sample.int(length(I_nodes), 1)
          i <- I_nodes[k]
          stopifnot(epi[i] == I)

          # remove from I (swap-delete)
          # remove from I_nodes
          last <- I_nodes[length(I_nodes)]

          I_nodes[k] <- last
          pos_in_I[last] <- k

          I_nodes <- I_nodes[-length(I_nodes)]
          pos_in_I[i] <- 0

          epi[i] <- R

          I_count <- I_count - 1
          R_count <- R_count + 1

          if (opinions[i] < 0) {

            I_red <- I_red - 1
            R_red <- R_red + 1

          } else {

            I_blue <- I_blue - 1
            R_blue <- R_blue + 1
          }

          S_red  <- total_red  - I_red - R_red
          S_blue <- total_blue - I_blue - R_blue

        }

      }
      
    # ---- record step ----
      event_counter <- event_counter + 1
      
      #=========================
      # Time
      #=========================
      time_hist <- c(time_hist, t)
      
      #=========================
      # Voter fractions
      #=========================
      if (event_counter %% record_every == 0) {
        frac_temp <- sapply(levels, function(op) {
          sum(opinions == op) / n
        })
        
        frac_mat <- cbind(frac_mat, frac_temp)
      }
      

      #=========================
      # Record epidemics
      #=========================
      if (isTRUE(input$runEpidemic)) {
      #--------------------------
      # Overall
      #--------------------------
      S_hist <- c(S_hist, sum(epi==S))
      I_hist <- c(I_hist, sum(epi==I))
      R_hist <- c(R_hist, sum(epi==R))

      #--------------------------
      # SIR by opinion camp
      #--------------------------
      camp <- ifelse(opinions < 0, "red", "blue")

      S_hist_red <- c(S_hist_red,  sum(epi[camp=="red"]  == S) )
      I_hist_red  <- c(I_hist_red, sum(epi[camp=="red"]  == I) )
      R_hist_red  <- c(R_hist_red, sum(epi[camp=="red"]  == R) )

      S_hist_blue <- c(S_hist_blue, sum(epi[camp=="blue"] == S) )
      I_hist_blue <- c(I_hist_blue, sum(epi[camp=="blue"] == I) )
      R_hist_blue <- c(R_hist_blue, sum(epi[camp=="blue"] == R) )

      #--------------------------
      # SIR by membership
      #--------------------------
      all_group_members <- unique(unlist(members))
      S_hist_grp <- c(S_hist_grp, sum(epi[all_group_members] == S))
      I_hist_grp <- c(I_hist_grp, sum(epi[all_group_members] == I))
      R_hist_grp <- c(R_hist_grp, sum(epi[all_group_members] == R))

      outsiders_all <- setdiff(1:length(epi), all_group_members)
      S_hist_out <- c(S_hist_out, sum(epi[outsiders_all] == S))
      I_hist_out <- c(I_hist_out, sum(epi[outsiders_all] == I))
      R_hist_out <- c(R_hist_out, sum(epi[outsiders_all] == R))
      
      }
      
# end gillespie
      
    }
    
    #=========================
    # Opinions
    #=========================
    opinion_history <- cbind(opinion_history, opinions)
    
    if (rig_dirty && isTRUE(input$show_rig)) {
      bipartite <- reconstruct_bipartite(members, n, m)
      RIG <- bipartite_to_rig(bipartite)
    }
    
  
    end_time <- Sys.time()
    comp_time <- end_time - start_time
    
    frac_mat <- t(frac_mat)
    colnames(frac_mat) <- levels
    frac_mat <- cbind(time_hist, frac_mat)
    
    # Final epidemic history
    SIR_df <- data.frame(time_hist,
                         S_hist/n,
                         I_hist/n,
                         R_hist/n)
    
    SIR_df_opinion_red <- data.frame(time_hist,
                                     S_hist_red/n,
                                     I_hist_red/n,
                                     R_hist_red/n)
    
    SIR_df_opinion_blue <- data.frame(time_hist,
                                      S_hist_blue/n,
                                      I_hist_blue/n,
                                      R_hist_blue/n)
    
    SIR_df_grp <- data.frame(time_hist,
                             S_hist_grp/n,
                             I_hist_grp/n,
                             R_hist_grp/n)
    
    SIR_df_out <- data.frame(time_hist,
                             S_hist_out/n,
                             I_hist_out/n,
                             R_hist_out/n)
    
    all(S_hist + I_hist + R_hist == n)
    
    # Output
    list(
      B0 = B0, B = bipartite,
      RIG = RIG, RIG0 = RIG0,
      opinions = opinions, #opinions0=opinions0,

      members0 = members0, frac_mat = frac_mat,
      opinion_history = opinion_history,
      num_opinions = num_opinions, members = members,

      SIR_df = SIR_df ,SIR_df_opinion_red = SIR_df_opinion_red,
      SIR_df_opinion_blue = SIR_df_opinion_blue,
      SIR_df_grp=SIR_df_grp, SIR_df_out = SIR_df_out,
      inf_time = inf_time, inf_camp = inf_camp,

      stop_reason = stop_reason, final_time = t, comp_time = comp_time
    )
  })

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


  #---------------
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

  output$bipartite0Plot <- renderPlot({     req(simData())           # simulation must exist
    req(input$show_rig0)     # checkbox must be TRUE;
    visual_bipartite(simData()$B0, simData()$opinion_history[,1], simData()$num_opinions) })

  output$bipartitePlot <- renderPlot({
    req(simData())           # simulation must exist
    req(input$show_rig)     # checkbox must be TRUE;
    visual_bipartite(simData()$B, simData()$opinion_history[,ncol( simData()$opinion_history)], simData()$num_opinions) })

  # Voter's dynamics
  output$fracPlot <- renderPlot({
    req(simData())

    visual_step_time(
      simData()$frac_mat,
      simData()$num_opinions
    )
  })

  output$histo <- renderPlot({ req(simData()); visual_histo(simData()$opinion_history, simData()$num_opinions) })
  output$histogramGroup0 <- renderPlot({ req(simData()); visual_histo_pergroup(simData()$opinion_history[,1], simData()$num_opinions, simData()$members0) })
  output$histogramGroup <- renderPlot({ req(simData()); visual_histo_pergroup(simData()$opinion_history[,ncol( simData()$opinion_history)], simData()$num_opinions, simData()$members) })

  # Epidemic layer
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

  # output$SIR6 <- renderPlot({
  #   req(simData())
  #   req(input$runEpidemic)
  #   visual_SIR(inf_time = simData()$inf_time, inf_camp = simData()$inf_camp)  })

  output$downloadAllPlots <- downloadHandler(

    filename = function() {
      paste0("all_plots_", Sys.Date(), ".zip")
    },

    content = function(file) {

      req(simData())

      # Temporary folder
      tmpdir <- tempdir()

      #==========================
      # Helper function
      #==========================
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

