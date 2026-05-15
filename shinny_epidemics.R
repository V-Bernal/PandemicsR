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
  titlePanel(h3("Voter - Schelling Multi-Membership with Epidemics Simulation")),


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
      checkboxInput("runSchelling", "Schelling model", value = TRUE),

      sliderInput("c_param", "Schelling: Edge addition rate param c", min = 0, step = 0.01, max = 1, value = 0.4),
      sliderInput("beta_plus", "Schelling: beta_plus", min = 0, step = 0.01, max = 1, value = 0.5),
      sliderInput("beta_minus", "Schelling: beta_minus", min = 0, step = 0.01, max = 1, value = 0.2),
      sliderInput("T_threshold", "Schelling: T_threshold", min = 0, step = 0.01, max = 1, value = 0.3)
      ),
      
      # Section 3: Voter
      wellPanel(
      h4("Voter's model"),
      checkboxInput("runVoter", "Voter model", value = TRUE),
      #sliderInput("kappa", "Voter: Poisson rate for opinion update kappa", min = 0, step = 0.01, max = 1, value = 0.3),
      sliderInput("gamma", "Voter: gamma", min = 0, step = 0.01, max = 100, value = 5),
      sliderInput("Numopinions", "Number of Opinions", min = 2, step = 2, max = 4, value = 2),

      # Section 4: Extremes
      #h5("Extremes"),
      sliderInput("alpha", "radicalization rate", min = 0, step = 0.01, max = 1, value = 0),
      sliderInput("alpha_deradicalization", "deradicalization rate*", min = 0, step = 0.01, max = 1, value = 0),
      h6("*Stubbornnes: set deradicalization to zero")
      ),

      # Section 5: Epidemics
      wellPanel(
      h4("Epidemics model"),
      checkboxInput("runEpidemic", "Epidemic model", value = TRUE),
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
      actionButton("runSim", "Run Simulation")
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
    # New Epidemic states
    #==========================
    # Each ind can be in one state S or I or R
    S <- 0; I <- 1; R <- 2

    # Initial Parameters (NEW)
    I0 <- input$I0 # number of initial infected
    beta_red  <- input$beta_red
    beta_blue <- input$beta_blue # infection rate blue opinion
    gamma_epi <- input$gamma_epi #recovery rate

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
      #group
      ids <- which(bipartite[, g] == 1)
      members[[g]] <- ids
      #opinion in group
      red_members[[g]] <- ids[opinions[ids] < 0] 
      blue_members[[g]] <- ids[opinions[ids] > 0 ]
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
    
    #=====================
    # Initialize epidemic (O(1) structure)
    #=====================
    
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
    
    # Global counters (O(1))
    I_count <- length(I_nodes)
    
    # # Camp-based infected counters
    # I_red  <- sum(opinions[I_nodes] < 0)
    # I_blue <- sum(opinions[I_nodes] > 0)
    # 
    # R_red <- 0
    # R_blue <- 0
    #=====================
    # Initialize epidemic history trackers
    #=====================
    
    # Initial time
    time_hist <- c(0)
    
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
    
    all_group_members <- unique(unlist(members))
    
    S_hist_grp <- c(sum(epi[all_group_members] == S))
    I_hist_grp <- c(sum(epi[all_group_members] == I))
    R_hist_grp <- c(sum(epi[all_group_members] == R))
    
    # Outsiders
    outsiders_all <- setdiff(seq_len(n), all_group_members)
    
    S_hist_out <- c(sum(epi[outsiders_all] == S))
    I_hist_out <- c(sum(epi[outsiders_all] == I))
    R_hist_out <- c(sum(epi[outsiders_all] == R))
    # #=====================
    # # Initialize epidemic
    # epi <- rep(S, n)
    # epi[sample(1:n, I0)] <- I
    # 
    # # Tracking Epidemics
    # time_hist <- c(0)
    # S_hist <- c(sum(epi==S))
    # I_hist <- c(sum(epi==I))
    # R_hist <- c(sum(epi==R))
    # 
    # camp <- ifelse(opinions < 0, "red", "blue")
    # S_hist_red <- c(sum(epi[camp=="red"]  == S))
    # I_hist_red  <- c(sum(epi[camp=="red"]  == I))
    # R_hist_red  <- c(sum(epi[camp=="red"]  == R))
    # 
    # S_hist_blue <- c(sum(epi[camp=="blue"] == S))
    # I_hist_blue  <- c(sum(epi[camp=="blue"] == I))
    # R_hist_blue  <- c(sum(epi[camp=="blue"] == R))
    # 
    # # SIR by group membership (aggregate)
    # # Flatten all members (unique individuals in any group)
    # all_group_members <- unique(unlist(members))
    # 
    # S_hist_grp <- sum(epi[all_group_members] == S)
    # I_hist_grp <- sum(epi[all_group_members] == I)
    # R_hist_grp <- sum(epi[all_group_members] == R)
    # 
    # # Outsiders
    # outsiders_all <- setdiff(1:length(epi), all_group_members)
    # 
    # S_hist_out <- sum(epi[outsiders_all] == S)
    # I_hist_out <- sum(epi[outsiders_all] == I)
    # R_hist_out <- sum(epi[outsiders_all] == R)

    #--time of infection----------
    inf_time <- c()
    inf_camp <- c()
    
    # # ---- ADD THIS ----
    # time_hist <- c()
    # S_hist <- c()
    # I_hist <- c()
    # R_hist <- c()
    # 
    # S_hist_red <- c(); I_hist_red <- c(); R_hist_red <- c()
    # S_hist_blue <- c(); I_hist_blue <- c(); R_hist_blue <- c()
    # S_hist_grp <- c(); I_hist_grp <- c(); R_hist_grp <- c()
    # S_hist_out <- c(); I_hist_out <- c(); R_hist_out <- c()
    # stopping time
    stop_reason <- "Reached maximum time"
    

    #==========================
    # Combined Voter and Schelling model dynamics Gillespie algorithm
    #==========================
    t <- 0
    while (t < t_max) {
      
      # stopping time
      if (all(epi == R)) {
        stop_reason <- "All individuals recovered (epidemic ended)"
        break
      }
      
      if (num_opinions == 4 && all(opinions %in% c(-2, 2))) {
        stop_reason <- "Full polarization (all extreme opinions)"
        break
      }
      
      
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

      #==========================
      # Epidemic rates (NEW)
      #==========================
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

      #==========================
      # New: extremist recruitment rates (and num of extremes)
      #==========================
      rate_extreme_red  <- numeric(m)
      rate_extreme_blue <- numeric(m)
      rate_deradicalize_red <- numeric(m)
      rate_deradicalize_blue <- numeric(m)


      for (g in seq_len(m)) {
        num_extreme_red   <- sum(opinions[members[[g]]] == -2)
        num_extreme_blue  <- sum(opinions[members[[g]]] ==  2)
        num_moderate_red  <- sum(opinions[members[[g]]] == -1)
        num_moderate_blue <- sum(opinions[members[[g]]] ==  1)
        
        Tot_g <- length(members[[g]])
        factor <- if (Tot_g > 0) 1 / Tot_g else 0 # avoid NAN or division by zero
        
        rate_extreme_red[g]  <- alpha * num_extreme_red  * num_moderate_red  * factor
        rate_extreme_blue[g] <- alpha * num_extreme_blue * num_moderate_blue * factor
        
        rate_deradicalize_red[g]  <- alpha_deradicalization * num_extreme_red  * num_moderate_red  * factor
        rate_deradicalize_blue[g] <- alpha_deradicalization * num_extreme_blue * num_moderate_blue * factor
        
      }
      #==========================

      # The rate at which the process leaves the state
      lambda_i <- voter_term + join_term + leaveR_rate + leaveB_rate +
        rate_extreme_red + rate_extreme_blue +
        rate_deradicalize_red + rate_deradicalize_blue

      lambda_tot <- sum(lambda_i) #+ infection_rate_tot + recovery_rate_tot
      if (lambda_tot <= 0) break


      #====================
      # Group selection at random
      #====================
      group_i <- if(m > 1) sample(1:m, 1L, prob = lambda_i / lambda_tot) else 1

      Ri_g <- Ri[group_i]; Bi_g <- Bi[group_i]
      voter_g <- voter_term[group_i]; join_g <- join_term[group_i]
      leaveR_g <- leaveR_rate[group_i]; leaveB_g <- leaveB_rate[group_i]

      rates_vec <- c(ifelse(voter_g>0, voter_g/2, 0),
                     ifelse(voter_g>0, voter_g/2, 0),
                     join_g, leaveR_g, leaveB_g)

      rates_vec <- c(rates_vec,
                     rate_extreme_red[group_i], # move 6
                     rate_extreme_blue[group_i] # move 7
      )

      rates_vec <- c(rates_vec,
                     rate_deradicalize_red[group_i], # move 8
                     rate_deradicalize_blue[group_i] # move 9
      )

      if (sum(rates_vec)<=0) next

      #====================
      # epidemic event selection
      #====================
      u_event <- runif(1)

      social_rate <- sum(lambda_i)
      epi_rate <- infection_rate_tot + recovery_rate_tot
      
      
      #================
      # Apply a random move in Gillespie time
      # Gillespie time
      #================
      dt <- rexp(1, lambda_tot)
      t <- t + dt
      if (t >= t_max) break

      if (u_event < social_rate / lambda_tot) {

      #================
      # Move step
      #================
      # 1: Red to Blue. (voter)
      # 2: Blue to Red. (voter)
      # 3: Add external to group
      # 4: Remove internal Red
      # 5: Remove internal Blue
      #============
      # 6 and 7: moderate turns extreme of same color (voter)
      # 8 and 9:  extreme to moderate of same color (voter)
      #============

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
      
      # Filter rates
      rates_sub <- rates_vec[enabled_moves]
      
      # Safety checks
      if (length(rates_sub) == 0 || sum(rates_sub) <= 0) next
      
      move <- sample(enabled_moves, 1L, prob = rates_sub / sum(rates_sub))
      
      #============

      #move <- sample(1:length(rates_vec), 1L, prob = rates_vec / sum(rates_vec))

      if (move==1 && Ri_g>0) { # Red to Blue
        chosen <- sample(red_members[[group_i]],1)
        opinions[chosen] <- +1
        
        if (epi[chosen] == I) {
          I_red  <- I_red - 1
          I_blue <- I_blue + 1
        }
        
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
        
        if (epi[chosen] == I) {
          I_blue <- I_blue - 1
          I_red  <- I_red + 1
        }
        
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
        if (opinions[chosen]< 0) {
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
      } else if (move == 6 &&
                 rate_extreme_red[group_i] > 0 &&
                 Ri[group_i] > 0) {

        moderates <- members[[group_i]][opinions[members[[group_i]]] == -1]
        extremes <- members[[group_i]][opinions[members[[group_i]]] == -2]

        if (length(moderates) > 0 && length(extremes) > 0) { #check if we need length(extremes) > 0
          chosen_moderate <- sample(moderates, 1)
          opinions[chosen_moderate] <- -2
          
          # The individual stays in the red camp, so red_members and Ri do not change.
        }
      } else if (move == 7 &&
                 rate_extreme_blue[group_i] > 0 &&
                 Bi[group_i] > 0) {

        moderates <- members[[group_i]][opinions[members[[group_i]]] == 1]
        extremes <- members[[group_i]][opinions[members[[group_i]]] == 2]
        if (length(moderates) > 0 && length(extremes) > 0) { #check if we need length(extremes) > 0
          chosen_moderate <- sample(moderates, 1)
          opinions[chosen_moderate] <- 2

          # The individual stays in the blue camp, so blue_members and Bi do not change.
        }
      } else if (move == 8 &&
                 rate_deradicalize_red[group_i] > 0 &&
                 Ri[group_i] > 0) {

        moderates <- members[[group_i]][opinions[members[[group_i]]] == -1]
        extremes <- members[[group_i]][opinions[members[[group_i]]] == -2]

        if (length(moderates) > 0 && length(extremes) > 0) {

          chosen_extreme <- sample(extremes, 1)
          opinions[chosen_extreme] <- -1
          
          # The individual stays in the red camp, so red_members and Ri do not change.
        }
      } else if (move == 9 &&
                 rate_deradicalize_blue[group_i] > 0 &&
                 Bi[group_i] > 0) {

        moderates <- members[[group_i]][opinions[members[[group_i]]] == 1]
        extremes <- members[[group_i]][opinions[members[[group_i]]] == 2]

        if (length(moderates) > 0 && length(extremes) > 0) {

          chosen_extreme <- sample(extremes, 1)
          opinions[chosen_extreme] <- 1
          
          # The individual stays in the blue camp, so blue_members and Bi do not change.
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
          # # probs <- infection_rates / sum(infection_rates)
          # # i <- sample(1:n, 1, prob = probs)
          # 
          # #=======================
          # inf_alias <- #build_alias(infection_rates / sum(infection_rates))
          # rec_alias <- # build_alias(recovery_rates / sum(recovery_rates))
          # 
          # if (runif(1) < infection_rate_tot / (infection_rate_tot + recovery_rate_tot)) {
          #   i <- alias_sample(inf_alias)
          # } else {
          #   i <- alias_sample(rec_alias)
          # }
          # #========================
          # 
          # # opinion at time of infection
          # camp_i <- ifelse(opinions[i] < 0, "red", "blue")
          # inf_time <- c(inf_time, t)
          # inf_camp <- c(inf_camp, camp_i)
          # 
          # epi[i] <- I
          
          # pick random susceptible
          if (length(S_nodes) == 0) next
          weights <- beta_vec[S_nodes]
          k <- sample.int(length(S_nodes), 1, prob = weights)
          i <- S_nodes[k]
          
          # opinion at time of infection
          camp_i <- ifelse(opinions[i] < 0, "red", "blue")
          inf_time <- c(inf_time, t)
          inf_camp <- c(inf_camp, camp_i)
          
          # remove from S (swap-delete)
          last <- S_nodes[length(S_nodes)]
          S_nodes[k] <- last
          pos_in_S[last] <- k
          S_nodes <- S_nodes[-length(S_nodes)]
          
          # add to I
          I_nodes <- c(I_nodes, i)
          pos_in_I[i] <- length(I_nodes)
          
          # update state
          epi[i] <- I
          
          # update counters
          S_count <- S_count - 1
          I_count <- I_count + 1
          
          if (opinions[i] < 0) I_red <- I_red + 1 else I_blue <- I_blue + 1

        } else {
          
          # -------- Recovery event --------
          # probs <- recovery_rates / sum(recovery_rates)
          # i <- sample(1:n, 1, prob = probs)
          # 
          # epi[i] <- R
          # pick random infected
          if (length(I_nodes) == 0) next
          k <- sample.int(length(I_nodes), 1)
          i <- I_nodes[k]
          
          # # opinion at time of infection
          # camp_i <- ifelse(opinions[i] < 0, "red", "blue")
          # inf_time <- c(inf_time, t)
          # inf_camp <- c(inf_camp, camp_i)
          
          # remove from I (swap-delete)
          last <- I_nodes[length(I_nodes)]
          I_nodes[k] <- last
          pos_in_I[last] <- k
          I_nodes <- I_nodes[-length(I_nodes)]
          
          # update state
          epi[i] <- R
          
          # update counters
          I_count <- I_count - 1
          R_count <- R_count + 1
          
          if (opinions[i] < 0) I_red <- I_red - 1 else I_blue <- I_blue - 1
          if (opinions[i] < 0) R_red <- R_red + 1 else R_blue <- R_blue + 1
          
          S_red  <- sum(opinions < 0) - I_red - R_red
          S_blue <- sum(opinions > 0) - I_blue - R_blue
          
        }
        
      }
      # ---- record step ----
      event_counter <- event_counter + 1

      #if (event_counter %% record_every == 0) {
      #  opinion_history <- cbind(opinion_history, opinions)
      #}

      # if (event_counter %% record_every == 0) {
      #   frac_temp <- sapply(levels, function(op) {
      #     sum(opinions == op)/n
      #   })
      #   frac_mat <- cbind(frac_mat, frac_temp)
      # }
      
      if (event_counter %% record_every == 0) {
        
        #=========================
        # Time
        #=========================
        time_hist <- c(time_hist, t)
        
        #=========================
        # Voter fractions
        #=========================
        frac_temp <- sapply(levels, function(op) {
          sum(opinions == op) / n
        })
        frac_mat <- cbind(frac_mat, frac_temp)
        
        #=========================
        # Epidemic: Overall
        #=========================
        S_hist <- c(S_hist, S_count)
        I_hist <- c(I_hist, I_count)
        R_hist <- c(R_hist, R_count)
        
        #=========================
        # Epidemic: Opinion camp (O(1))
        #=========================
        S_red  <- sum(opinions < 0) - I_red - R_red
        S_blue <- sum(opinions > 0) - I_blue - R_blue
        
        S_hist_red  <- c(S_hist_red,  S_red)
        I_hist_red  <- c(I_hist_red,  I_red)
        R_hist_red  <- c(R_hist_red,  R_red)
        
        S_hist_blue <- c(S_hist_blue, S_blue)
        I_hist_blue <- c(I_hist_blue, I_blue)
        R_hist_blue <- c(R_hist_blue, R_blue)
        
        #=========================
        # Epidemic: Group membership
        #=========================
        all_group_members <- unique(unlist(members))
        
        S_hist_grp <- c(S_hist_grp, sum(epi[all_group_members] == S))
        I_hist_grp <- c(I_hist_grp, sum(epi[all_group_members] == I))
        R_hist_grp <- c(R_hist_grp, sum(epi[all_group_members] == R))
        
        outsiders_all <- setdiff(seq_len(n), all_group_members)
        
        S_hist_out <- c(S_hist_out, sum(epi[outsiders_all] == S))
        I_hist_out <- c(I_hist_out, sum(epi[outsiders_all] == I))
        R_hist_out <- c(R_hist_out, sum(epi[outsiders_all] == R))
      }


    }

    # Final opinion history
    opinion_history <- cbind(opinion_history, opinions)
    frac_mat <- t(frac_mat)
    colnames(frac_mat) <- levels
    frac_mat <-cbind(time_hist, frac_mat)

    if (rig_dirty && isTRUE(input$show_rig)) {
      bipartite <- reconstruct_bipartite(members, n, m)
      RIG <- bipartite_to_rig(bipartite)
    }

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
    
    end_time <- Sys.time()
    comp_time <- end_time - start_time

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
           round(simData()$comp_time, 3),
           " seconds")
  })
  
  
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

  # Voter's dynamics
  output$fracPlot <- renderPlot({
    req(simData())
    
    visual_step_time(
      simData()$frac_mat,
      simData()$num_opinions
    )
  })  
  
  output$histo <- renderPlot({ req(simData()); visual_histo(simData()$opinion_history, simData()$num_opinions) })
  #output$heatmapPlot <- renderPlot({ req(simData()); heatmapPlot(simData()$opinion_history, simData()$num_opinions) })
  #output$histogramGroup0 <- renderPlot({ req(simData()); visual_histo_pergroup(simData()$opinion_history[,1], simData()$num_opinions, simData()$B0) })
  #output$histogramGroup <- renderPlot({ req(simData()); visual_histo_pergroup(simData()$opinion_history[,ncol( simData()$opinion_history)], simData()$num_opinions, simData()$B) })
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

  }

#==========================
# --- Run App ---
#==========================
shinyApp(ui, server)
