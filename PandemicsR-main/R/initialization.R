#==========================
# Initialize Network
#==========================
init_network <- function(params){
  # Homogeneous weights
  ind_w <- rep(params$lambda, params$n)
  grp_w <- rep(params$lambda * params$n / params$m, params$m)
  # bipartite and RIG
  B0 <- generate_bipartite(
    params$n,
    params$m,
    ind_w,
    grp_w
  )

  list(
    B0 = B0,
    RIG0 = bipartite_to_rig(B0)
  )
}

#==========================
# Initialize Opinion trackers
#==========================
# Voter tracker
init_opinions <- function(params){

  opinions <- initialize_opinions_multi(
    params$n,
    params$num_opinions
  )

  opinion_history <- data.frame(
    matrix(
      opinions,
      nrow = params$n,
      ncol = 2
    )
  )

  list(
    opinions = opinions,
    opinion_history = opinion_history
  )
}

#==========================
# Initialize Groups
#==========================
init_groups <- function(B0, opinions, params){

  m <- params$m
  n <- params$n

  members <- vector("list", m)
  red_members <- vector("list", m)
  blue_members <- vector("list", m)
  outsiders <- vector("list", m)

  groups_of_individual <- vector("list", n)

  # Group membership trackers
  for (g in seq_len(m)) {

    ids <- which(B0[, g] == 1)

    members[[g]] <- ids

    red_members[[g]] <- ids[opinions[ids] == -1]

    blue_members[[g]] <- ids[opinions[ids] == 1]

    outsiders[[g]] <- setdiff(seq_len(n), ids)

    for (i in ids) {
      groups_of_individual[[i]] <-
        c(groups_of_individual[[i]], g)
    }
  }


  # Moderate counts
  Ri <- sapply(red_members, length)
  Bi <- sapply(blue_members, length)


  # Extreme counts
  Ei_red <- integer(m)
  Ei_blue <- integer(m)

  for (g in seq_len(m)) {

    ids <- members[[g]]

    Ei_red[g] <-
      sum(opinions[ids] == -2)

    Ei_blue[g] <-
      sum(opinions[ids] == 2)
  }


  list(
    members = members,
    red_members = red_members,
    blue_members = blue_members,
    outsiders = outsiders,
    groups_of_individual = groups_of_individual,
    Ri = Ri,
    Bi = Bi,
    Ei_red = Ei_red,
    Ei_blue = Ei_blue
  )
}

#==========================
# Initialize Epidemics
#==========================
init_epidemic <- function(params, opinions){

  if (!params$runEpidemic) {
    return(NULL)
  }

  n <- params$n

  # States
  S <- 0
  I <- 1
  R <- 2

  # Initial infected
  I0 <- floor(params$I0 * n)

  epi <- rep(S, n)

  initial_infected <- sample.int(n, I0)

  epi[initial_infected] <- I


  # Explicit sets
  S_nodes <- setdiff(seq_len(n), initial_infected)
  I_nodes <- initial_infected
  R_nodes <- integer(0)

  S_count <- length(S_nodes)
  I_count <- length(I_nodes)
  R_count <- 0

  # Maintain susceptible sets by opinion
  S_red_nodes  <- S_nodes[opinions[S_nodes] < 0]
  S_blue_nodes <- S_nodes[opinions[S_nodes] > 0]


  # Positions for O(1) removal
  pos_in_S <- integer(n)
  pos_in_I <- integer(n)
  pos_in_R <- integer(n)

  for(k in seq_along(S_nodes)){
    pos_in_S[S_nodes[k]] <- k
  }

  for(k in seq_along(I_nodes)){
    pos_in_I[I_nodes[k]] <- k
  }


  # Opinion camp counters
  total_red  <- sum(opinions < 0)
  total_blue <- sum(opinions > 0)

  I_red  <- sum(opinions[I_nodes] < 0)
  I_blue <- sum(opinions[I_nodes] > 0)

  R_red  <- 0
  R_blue <- 0

  S_red  <- total_red  - I_red  - R_red
  S_blue <- total_blue - I_blue - R_blue


  # infection history
  inf_time <- numeric(0)
  inf_camp <- character(0)


  list(

    S = S,
    I = I,
    R = R,

    epi = epi,

    S_nodes = S_nodes,
    I_nodes = I_nodes,
    R_nodes = R_nodes,

    S_red_nodes = S_red_nodes,
    S_blue_nodes = S_blue_nodes,

    S_count = S_count,
    I_count = I_count,
    R_count = R_count,

    pos_in_S = pos_in_S,
    pos_in_I = pos_in_I,
    pos_in_R = pos_in_R,

    total_red = total_red,
    total_blue = total_blue,

    S_red = S_red,
    I_red = I_red,
    R_red = R_red,

    S_blue = S_blue,
    I_blue = I_blue,
    R_blue = R_blue,

    inf_time = inf_time,
    inf_camp = inf_camp
  )
}


#==========================
# Format Output
#==========================
finalize_output <- function(
    state,
    trackers,
    params,
    network_state,
    stop_reason,
    t,
    comp_time,
    event_counter){

  # Reconstruct RIG if needed
  B <- network_state$B0
  RIG <- network_state$RIG0

  if (isTRUE(state$rig_dirty)) {

    B <- reconstruct_bipartite(
      state$members,
      params$n,
      params$m
    )

    RIG <- bipartite_to_rig(B)
  }


  tracker_output <- finalize_trackers(
    trackers,
    params
  )


  #=========================
  # Return final object
  #=========================

  list(

    #-----------------------
    # Network
    #-----------------------
    B0 = network_state$B0,
    B  = B,

    RIG0 = network_state$RIG0,
    RIG  = RIG,


    #-----------------------
    # Final opinion state
    #-----------------------
    num_opinions = params$num_opinions,
    opinions = state$opinions,
    opinion_history = trackers$opinion_history,

    members0 = state$members0,
    members  = state$members,


    #-----------------------
    # Opinion history
    #-----------------------
    frac_mat =
      tracker_output$frac_mat,


    #-----------------------
    # Epidemic histories
    #-----------------------
    SIR_df =
      tracker_output$SIR_df,

    SIR_df_opinion_red =
      tracker_output$SIR_df_opinion_red,

    SIR_df_opinion_blue =
      tracker_output$SIR_df_opinion_blue,

    SIR_df_grp =
      tracker_output$SIR_df_grp,

    SIR_df_out =
      tracker_output$SIR_df_out,


    #-----------------------
    # Infection events
    #-----------------------
    inf_time =
      state$inf_time,

    inf_camp =
      state$inf_camp,


    #-----------------------
    # Metadata
    #-----------------------
    stop_reason = stop_reason,

    final_time = t,

    comp_time = comp_time,

    event_counter = event_counter
  )

}
