apply_epidemic_event <- function(state, rates, params, t) {

  #==========================
  # Safety check
  #==========================

  if (rates$epi_rate <= 0)
    return(state)


  #==========================
  # Infection vs recovery
  #==========================

  if (stats::runif(1) < rates$infection_rate_tot / rates$epi_rate) {

    # Infection event
    if (length(state$S_nodes) == 0)
      return(state)

    weights <- rates$beta_vec[state$S_nodes]

    if (sum(weights) <= 0)
      return(state)


    # choose susceptible
    k <- sample.int(
      length(state$S_nodes),
      1,
      prob = weights
    )

    i <- state$S_nodes[k]

    #stopifnot(state$epi[i] == state$S)


    #--------------------------
    # infection history
    #--------------------------

    camp_i <- ifelse(
      state$opinions[i] < 0,
      "red",
      "blue"
    )

    state$inf_time <- c(
      state$inf_time,
      t
    )

    state$inf_camp <- c(
      state$inf_camp,
      camp_i
    )


    #--------------------------
    # remove from S_nodes
    # swap-delete
    #--------------------------

    last <- state$S_nodes[length(state$S_nodes)]

    state$S_nodes[k] <- last
    state$pos_in_S[last] <- k

    state$S_nodes <- state$S_nodes[-length(state$S_nodes)]

    state$pos_in_S[i] <- 0


    #--------------------------
    # add to I_nodes
    #--------------------------

    state$I_nodes <- c(
      state$I_nodes,
      i
    )

    state$pos_in_I[i] <- length(state$I_nodes)


    #--------------------------
    # update states
    #--------------------------

    state$epi[i] <- state$I

    state$S_count <- state$S_count - 1
    state$I_count <- state$I_count + 1


    #--------------------------
    # update camps
    #--------------------------

    if (state$opinions[i] < 0) {

      state$I_red <- state$I_red + 1

    } else {

      state$I_blue <- state$I_blue + 1

    }


    state$S_red <-
      state$total_red -
      state$I_red -
      state$R_red

    state$S_blue <-
      state$total_blue -
      state$I_blue -
      state$R_blue


  } else {

    # Recovery event
    if (length(state$I_nodes) == 0)
      return(state)


    k <- sample.int(
      length(state$I_nodes),
      1
    )

    i <- state$I_nodes[k]

    #stopifnot(state$epi[i] == state$I)


    #--------------------------
    # remove from I_nodes
    #--------------------------

    last <- state$I_nodes[length(state$I_nodes)]

    state$I_nodes[k] <- last
    state$pos_in_I[last] <- k

    state$I_nodes <- state$I_nodes[-length(state$I_nodes)]

    state$pos_in_I[i] <- 0

    #--------------------------
    # add to R_nodes
    #--------------------------

    state$R_nodes <- c(state$R_nodes, i)
    state$pos_in_R[i] <- length(state$R_nodes)


    #--------------------------
    # update epidemic state
    #--------------------------

    state$epi[i] <- state$R

    state$I_count <- state$I_count - 1
    state$R_count <- state$R_count + 1


    #--------------------------
    # update camps
    #--------------------------

    if (state$opinions[i] < 0) {

      state$I_red <- state$I_red - 1
      state$R_red <- state$R_red + 1

    } else {

      state$I_blue <- state$I_blue - 1
      state$R_blue <- state$R_blue + 1

    }


    state$S_red <-
      state$total_red -
      state$I_red -
      state$R_red

    state$S_blue <-
      state$total_blue -
      state$I_blue -
      state$R_blue

  }


  return(state)
}
