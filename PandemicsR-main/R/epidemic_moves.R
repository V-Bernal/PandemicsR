apply_epidemic_event <- function(state, rates, params, t) {

  if (rates$epi_rate <= 0)
    return(state)

    #==========================
    # Select Infection or Recovery
    #==========================

    if (stats::runif(1) < rates$infection_rate_tot / rates$epi_rate) {

      #  1. Infection
      # Choose infection camp
      if (rates$infection_rate_tot <= 0)
        return(state)

      if (stats::runif(1) <
          rates$infection_red_rate / rates$infection_rate_tot) {

        if (length(state$S_red_nodes) == 0)
          return(state)

        # remove from red susceptible set
        i <- sample(length(state$S_red_nodes), 1)

        state$S_red_nodes <-
          state$S_red_nodes[state$S_red_nodes != i]

      } else {

        # remove from blue susceptible set
        if (length(state$S_blue_nodes) == 0)
          return(state)

        i <- sample(length(state$S_blue_nodes), 1)


        state$S_blue_nodes <-
          state$S_blue_nodes[state$S_blue_nodes != i]

      }

        # Infection history
        state$inf_time <- c(state$inf_time, t)
        state$inf_camp <- c(
          state$inf_camp,
          ifelse(state$opinions[i] < 0, "red", "blue")
      )

      # Remove from S_nodes (keep your existing swap-delete)
      k <- state$pos_in_S[i]

      last <- state$S_nodes[length(state$S_nodes)]

      state$S_nodes[k] <- last
      state$pos_in_S[last] <- k

      state$S_nodes <- state$S_nodes[-length(state$S_nodes)]
      state$pos_in_S[i] <- 0

      # Add to infected
      state$I_nodes <- c(state$I_nodes, i)
      state$pos_in_I[i] <- length(state$I_nodes)

      state$epi[i] <- state$I

      state$S_count <- state$S_count - 1
      state$I_count <- state$I_count + 1

      if (state$opinions[i] < 0) {

        state$I_red <- state$I_red + 1

      } else {

        state$I_blue <- state$I_blue + 1

      }

      state$S_red  <- state$total_red  - state$I_red  - state$R_red
      state$S_blue <- state$total_blue - state$I_blue - state$R_blue

    } else {

      # 2. Recovery
      # Choose recover
      if (length(state$I_nodes) == 0)
        return(state)

      k <- sample.int(
        length(state$I_nodes),
        1
      )

      i <- state$I_nodes[k]

      stopifnot(state$epi[i] == state$I)

      # remove from I_nodes
      last <- state$I_nodes[length(state$I_nodes)]

      state$I_nodes[k] <- last
      state$pos_in_I[last] <- k

      state$I_nodes <- state$I_nodes[-length(state$I_nodes)]

      state$pos_in_I[i] <- 0

      # add to R_nodes
      state$R_nodes <- c(state$R_nodes, i)
      state$pos_in_R[i] <- length(state$R_nodes)


      # update epidemic state
      state$epi[i] <- state$R

      state$I_count <- state$I_count - 1
      state$R_count <- state$R_count + 1

      # update camps
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
