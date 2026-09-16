#' epidemic_event
#'
#' @param params Simulation parameters.
#' @param opinions
#' @param B0
#' @param params
#' @param state,
#' @param trackers
#' @param params
#' @param network_state
#' @param stop_reason
#' @param t
#' @param comp_time
#' @param event_counter
#' @export
apply_epidemic_event <- function(state, rates, params, t) {

  if (rates$epi_rate <= 0)
    return(state)

    # Select Infection or Recovery

    if (stats::runif(1) < rates$infection_rate_tot / rates$epi_rate) {

      #  1. Infection

      #  Choose opinion
      if (rates$infection_rate_tot <= 0)
        return(state)

      if (stats::runif(1) <
          rates$infection_red_rate / rates$infection_rate_tot) {

        if (length(state$S_red_nodes) == 0)
          return(state)

        # remove i from red susceptibles
        k <- sample.int(length(state$S_red_nodes), 1)
        i <- state$S_red_nodes[k]

        # Remove from S_red using swap-and-pop
        last <- length(state$S_red_nodes)
        last_node <- state$S_red_nodes[last]

        state$S_red_nodes[k] <- last_node
        state$S_red_pos[last_node] <- k
        state$S_red_pos[i] <- 0L
        length(state$S_red_nodes) <- last - 1L

        # Add to I_red
        state$I_red_nodes <- c(state$I_red_nodes, i)
        state$I_red_pos[i] <- length(state$I_red_nodes)

        state$I_red <- state$I_red + 1
        state$S_red <- state$total_red - state$I_red - state$R_red
        # k <- sample.int(length(state$S_red_nodes), 1)
        # i <- state$S_red_nodes[k]
        #
        #   #--- new: swap and pop
        #   # state$S_red_nodes <-
        #   #  state$S_red_nodes[state$S_red_nodes != i]
        #   last <- length(state$S_red_nodes)
        #   state$S_red_nodes[k] <- state$S_red_nodes[last]
        #   state$S_red_nodes <- state$S_red_nodes[-last]
        #   #-----
        #
        # state$I_red <- state$I_red + 1
        # state$S_red  <- state$total_red  - state$I_red  - state$R_red
        # state$I_red_nodes <- c(state$I_red_nodes, i)

      } else {

        # 1.2.  remove from blue susceptible set
        if (length(state$S_blue_nodes) == 0)
          return(state)

        k <- sample.int(length(state$S_blue_nodes), 1)
        i <- state$S_blue_nodes[k]

        # Remove from S_blue using swap-and-pop
        last <- length(state$S_blue_nodes)
        last_node <- state$S_blue_nodes[last]

        state$S_blue_nodes[k] <- last_node
        state$S_blue_pos[last_node] <- k
        state$S_blue_pos[i] <- 0L
        length(state$S_blue_nodes) <- last - 1L

        # Add to I_blue
        state$I_blue_nodes <- c(state$I_blue_nodes, i)
        state$I_blue_pos[i] <- length(state$I_blue_nodes)

        state$I_blue <- state$I_blue + 1
        state$S_blue <- state$total_blue - state$I_blue - state$R_blue

        # k <- sample.int(length(state$S_blue_nodes), 1)
        # i <- state$S_blue_nodes[k]
        #
        # # new: swap and pop
        # # state$S_blue_nodes <-
        # #  state$S_blue_nodes[state$S_blue_nodes != i]
        # last <- length(state$S_blue_nodes)
        # state$S_blue_nodes[k] <- state$S_blue_nodes[last]
        # state$S_blue_nodes <- state$S_blue_nodes[-last]
        # #
        #
        # state$I_blue <- state$I_blue + 1
        # state$S_blue <- state$total_blue - state$I_blue - state$R_blue
        # state$I_blue_nodes <- c(state$I_blue_nodes, i)

      }

        # Infection history
        state$inf_time <- c(state$inf_time, t)
        state$inf_camp <- c(
          state$inf_camp,
          ifelse(state$opinions[i] < 0, "red", "blue"))

        state$S_count <- state$S_count - 1
        state$I_count <- state$I_count + 1

        state$I_nodes <- c(state$I_nodes, i)
        state$I_pos[i] <- length(state$I_nodes)

       # Remove from S_nodes (keep your existing swap-delete)
       # state$S_nodes <- state$S_nodes[state$S_nodes != i]
        idx <- state$S_pos[i]
        last <- length(state$S_nodes)
        last_node <- state$S_nodes[last]

        state$S_nodes[idx] <- last_node
        state$S_pos[last_node] <- idx
        state$S_pos[i] <- 0L
        length(state$S_nodes) <- last - 1L

      state$epi[i] <- state$I

    } else {

      # 2. Recovery
      # Choose recover
      if (length(state$I_nodes) == 0)
        return(state)

      k <- sample.int(length(state$I_nodes), 1)
      i <- state$I_nodes[k]

      stopifnot(state$epi[i] == state$I)

      state$I_count <- state$I_count - 1
      state$R_count <- state$R_count + 1

      # remove from I_nodes
      #state$I_nodes <- state$I_nodes[state$I_nodes != i]
      idx <- state$I_pos[i]
      last <- length(state$I_nodes)
      last_node <- state$I_nodes[last]

      state$I_nodes[idx] <- last_node
      state$I_pos[last_node] <- idx
      state$I_pos[i] <- 0L
      length(state$I_nodes) <- last - 1L

      # add to R_nodes
      #state$R_nodes <- c(state$R_nodes, i)
      state$R_nodes <- c(state$R_nodes, i)
      state$R_pos[i] <- length(state$R_nodes)

      # update epidemic state
      state$epi[i] <- state$R

      # update camps
      if (state$opinions[i] < 0) {

        state$I_red <- state$I_red - 1
        state$R_red <- state$R_red + 1

        # state$I_red_nodes <-
        # state$I_red_nodes[state$I_red_nodes != i]
        #
        # state$R_red_nodes <-
        #   c(state$R_red_nodes, i)
        # Remove from I_red using swap-and-pop
        idx <- state$I_red_pos[i]
        last <- length(state$I_red_nodes)
        last_node <- state$I_red_nodes[last]

        state$I_red_nodes[idx] <- last_node
        state$I_red_pos[last_node] <- idx
        state$I_red_pos[i] <- 0L
        length(state$I_red_nodes) <- last - 1L

        # Add to R_red
        state$R_red_nodes <- c(state$R_red_nodes, i)
        state$R_red_pos[i] <- length(state$R_red_nodes)

        # cat(
        #   "\nRECOVERY RED:",
        #   "i =", i,
        #   "R_red_nodes =", paste(state$R_red_nodes, collapse = ","),
        #   "R_red_pos[i] =", state$R_red_pos[i],
        #   "\n"
        # )

        stopifnot(
          i %in% state$R_red_nodes,
          state$R_red_pos[i] > 0L,
          state$R_red_nodes[state$R_red_pos[i]] == i
        )

      } else {

        state$I_blue <- state$I_blue - 1
        state$R_blue <- state$R_blue + 1

        # state$I_blue_nodes <- state$I_blue_nodes[state$I_blue_nodes != i]
        # state$R_blue_nodes <- c(state$R_blue_nodes, i)
        # Remove from I_blue using swap-and-pop
        idx <- state$I_blue_pos[i]
        last <- length(state$I_blue_nodes)
        last_node <- state$I_blue_nodes[last]

        state$I_blue_nodes[idx] <- last_node
        state$I_blue_pos[last_node] <- idx
        state$I_blue_pos[i] <- 0L
        length(state$I_blue_nodes) <- last - 1L

        # Add to R_blue
        state$R_blue_nodes <- c(state$R_blue_nodes, i)
        state$R_blue_pos[i] <- length(state$R_blue_nodes)

        # cat(
        #   "\nRECOVERY BLUE:",
        #   "i =", i,
        #   "R_blue_nodes =", paste(state$R_blue_nodes, collapse = ","),
        #   "R_blue_pos[i] =", state$R_blue_pos[i],
        #   "\n"
        # )

        stopifnot(
          i %in% state$R_blue_nodes,
          state$R_blue_pos[i] > 0L,
          state$R_blue_nodes[state$R_blue_pos[i]] == i
        )

      }


    }

  check_epidemic_position_indexes(state)
  check_epidemic_camp_membership(state)
  return(state)
}

