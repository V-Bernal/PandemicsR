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
      checkboxInput("runStubborn", "Stubborn", value = TRUE),
      sliderInput("alpha", "radicalization rate", min = 0, step = 0.01, max = 1, value = 0),
      sliderInput("alpha_deradicalization", "deradicalization rate", min = 0, step = 0.01, max = 1, value = 0),


      # Section 5: Epidemics
      h4("Epidemic's parameters"),
      checkboxInput("runEpidemic", "Epidemic model", value = TRUE),
      sliderInput("I0", "Number of Infected", min = 1, step = 1, max = 20, value = 0),
      sliderInput("gamma_epi", "Epidemic: recovery rate", min = 0, step = 0.01, max = 1, value = 0.3),
      sliderInput("beta_red", "Epidemic: infection rate red", min = 0, step = 0.01, max = 1, value = 0.5),
      sliderInput("beta_blue", "Epidemic: infection rate blue", min = 0, step = 0.01, max = 1, value = 0.2),

      # Section 6: Visualization
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

                 h4("Bipartite graph"),
                 plotOutput("bipartite0Plot", height = "500px"),
                 plotOutput("bipartitePlot", height = "500px"),
        ),

        tabPanel("Voter's dynamics",

                 h4("Time evolution of opinions"),
                 plotOutput("fracPlot", width = "500px", height = "400px"),
                 #plotOutput("heatmapPlot", height = "350px")

                 h4("Overall opinions"),
                 plotOutput("histo"),

                 h4("Group-wise initial opinions"),
                 plotOutput("histogramGroup0"),

                 h4("Group-wise final opinions"),
                 plotOutput("histogramGroup")
        ),
        tabPanel("Epidemic's dynamics",

                 h4("SIR overall"),
                 plotOutput("SIR1"),

                 h4("SIR red opinion"),
                 plotOutput("SIR2"),

                 h4("SIR blue opinion"),
                 plotOutput("SIR3"),

                 h4("SIR members"),
                 plotOutput("SIR4"),

                 h4("SIR isolated"),
                 plotOutput("SIR5")#,
#
#                  h4("SIR cummulative infections"),
#                  plotOutput("SIR6")
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

    #=====================
    # Initialize epidemic
    epi <- rep(S, n)
    epi[sample(1:n, I0)] <- I

    # Tracking Epidemics
    time_hist <- c(0)
    S_hist <- c(sum(epi==S))
    I_hist <- c(sum(epi==I))
    R_hist <- c(sum(epi==R))

    camp <- ifelse(opinions < 0, "red", "blue")
    S_hist_red <- c(sum(epi[camp=="red"]  == S))
    I_hist_red  <- c(sum(epi[camp=="red"]  == I))
    R_hist_red  <- c(sum(epi[camp=="red"]  == R))

    S_hist_blue <- c(sum(epi[camp=="blue"] == S))
    I_hist_blue  <- c(sum(epi[camp=="blue"] == I))
    R_hist_blue  <- c(sum(epi[camp=="blue"] == R))

    # SIR by group membership (aggregate)
    # Flatten all members (unique individuals in any group)
    all_group_members <- unique(unlist(members))

    S_hist_grp <- sum(epi[all_group_members] == S)
    I_hist_grp <- sum(epi[all_group_members] == I)
    R_hist_grp <- sum(epi[all_group_members] == R)

    # Outsiders
    outsiders_all <- setdiff(1:length(epi), all_group_members)

    S_hist_out <- sum(epi[outsiders_all] == S)
    I_hist_out <- sum(epi[outsiders_all] == I)
    R_hist_out <- sum(epi[outsiders_all] == R)

    #--time of infection----------
    inf_time <- c()
    inf_camp <- c()

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

      #==========================
      # Epidemic rates (NEW)
      #==========================
      I_count <- sum(epi == I)
      prevalence <- I_count / n # global = random mixing

      # Opinion camp-based infection rates
      beta_vec <- ifelse(opinions < 0, beta_red, beta_blue)

      # Infection: only for S individuals
      infection_rates <- ifelse(epi == S, beta_vec * prevalence, 0)

      # Recovery: only for I individuals
      recovery_rates <- ifelse(epi == I, gamma_epi, 0)

      # Aggregate epidemic rates
      infection_rate_tot <- sum(infection_rates)
      recovery_rate_tot  <- sum(recovery_rates)


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
        #Tot_g <- length(members[[g]])
        #valid <- Tot_g >= 2
        rate_extreme_red[g]  <- alpha * num_extreme_red  * num_moderate_red #/ Tot_g
        rate_extreme_blue[g] <- alpha * num_extreme_blue * num_moderate_blue #/ Tot_g

        # same
        rate_deradicalize_red[g] <- alpha_deradicalization * num_extreme_red  * num_moderate_red
        rate_deradicalize_blue[g] <- alpha_deradicalization * num_extreme_blue * num_moderate_blue
      }
      #==========================

      # The rate at which the process leaves the state
      lambda_i <- voter_term + join_term + leaveR_rate + leaveB_rate +
        rate_extreme_red + rate_extreme_blue +
        rate_deradicalize_red + rate_deradicalize_blue

      lambda_tot <- sum(lambda_i) + infection_rate_tot + recovery_rate_tot
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

      if (u_event < social_rate / lambda_tot) {


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
      # 6 and 7: moderate turns extreme of same color (voter)
      # 8 and 9:  extreme to moderate of same color (voter)
      #============

      #============
      # Trim rate vec according to the number of moves in simulation
      # Voter: moves 1 and 2
      # Schelling: moves 3 to 5
      trimmer <- c(rep(input$runVoter, 2), rep(input$runSchelling, 3), rep(num_opinions == 4, 2), rep(input$Stubborn != T, 2))
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
      } else if (move == 6 &&
                 rate_extreme_red[group_i] > 0 &&
                 Ri[group_i] > 0) {

        moderates <- members[[group_i]][opinions[members[[group_i]]] == -1]
        extremes <- members[[group_i]][opinions[members[[group_i]]] == -2]

        if (length(moderates) > 0 && length(extremes) > 0) { #check if we need length(extremes) > 0
          chosen_moderate <- sample(moderates, 1)
          opinions[chosen_moderate] <- -2

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
                 rate_extreme_blue[group_i] > 0 &&
                 Bi[group_i] > 0) {

        moderates <- members[[group_i]][opinions[members[[group_i]]] == 1]
        extremes <- members[[group_i]][opinions[members[[group_i]]] == 2]
        if (length(moderates) > 0 && length(extremes) > 0) { #check if we need length(extremes) > 0
          chosen_moderate <- sample(moderates, 1)
          opinions[chosen_moderate] <- 2

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
      } else if (move == 8 &&
                 rate_deradicalize_red[group_i] > 0 &&
                 Ri[group_i] > 0) {

        moderates <- members[[group_i]][opinions[members[[group_i]]] == -1]
        extremes <- members[[group_i]][opinions[members[[group_i]]] == -2]

        if (length(moderates) > 0 && length(extremes) > 0) {

          chosen_extreme <- sample(extremes, 1)
          opinions[chosen_extreme] <- -1

          # UPDATE ALL GROUPS OF THE INDIVIDUAL
          for (g in groups_of_individual[[chosen_extreme]]) {

            # add to red_members (moderates)
            red_members[[g]] <- c(red_members[[g]], chosen_extreme)

            # update count
            Ri[g] <- Ri[g] + 1
          }
        }
      } else if (move == 9 &&
                 rate_deradicalize_blue[group_i] > 0 &&
                 Bi[group_i] > 0) {

        moderates <- members[[group_i]][opinions[members[[group_i]]] == -1]
        extremes <- members[[group_i]][opinions[members[[group_i]]] == -2]

        if (length(moderates) > 0 && length(extremes) > 0) {

          chosen_extreme <- sample(extremes, 1)
          opinions[chosen_extreme] <- 1

          # UPDATE ALL GROUPS OF THE INDIVIDUAL
          for (g in groups_of_individual[[chosen_extreme]]) {

            # add to blue_members (moderates)
            blue_members[[g]] <- c(blue_members[[g]], chosen_extreme)

            # update count
            Bi[g] <- Bi[g] + 1
          }
        }
      } else {

        #==========================
        # Section 5. Epidemic event
        #==========================

        if ((infection_rate_tot + recovery_rate_tot) == 0) next

        # Decide infection vs recovery
        if (runif(1) < infection_rate_tot / (infection_rate_tot + recovery_rate_tot)) {

          # -------- Infection event --------
          probs <- infection_rates / sum(infection_rates)
          i <- sample(1:n, 1, prob = probs)

          # opinion at time of infection
          camp_i <- ifelse(opinions[i] < 0, "red", "blue")
          inf_time <- c(inf_time, t)
          inf_camp <- c(inf_camp, camp_i)

          epi[i] <- I

        } else {

          # -------- Recovery event --------
          probs <- recovery_rates / sum(recovery_rates)
          i <- sample(1:n, 1, prob = probs)

          epi[i] <- R
        }
      }

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

    #=========================
    # Record epidemics
    #=========================

    #--------------------------
    # Overall
    #--------------------------
    time_hist <- c(time_hist, t)
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

    # Final opinion history
    opinion_history <- cbind(opinion_history, opinions)
    frac_mat <- t(frac_mat)
    colnames(frac_mat) <- levels

    # if (rig_dirty) {
    #   bipartite <- reconstruct_bipartite(members, n, m)
    #   RIG <- bipartite_to_rig(bipartite)
    # }

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
      inf_time = inf_time, inf_camp = inf_camp
    )
  })

  #==========================
  # Section 6: Visualization
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

  # Voter's dynamics
  output$fracPlot <- renderPlot({ req(simData()); visual_step_time(simData()$frac_mat, simData()$num_opinions, length(simData()$opinion_history[,1])) })
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
