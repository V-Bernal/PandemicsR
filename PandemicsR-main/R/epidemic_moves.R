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
      #  1.1.Choose opinion
      if (rates$infection_rate_tot <= 0)
        return(state)

      if (stats::runif(1) <
          rates$infection_red_rate / rates$infection_rate_tot) {

        if (length(state$S_red_nodes) == 0)
          return(state)

        # 1.2. remove i from red susceptibles (swap-and-pop O(1))
        k <- sample.int(length(state$S_red_nodes), 1)
        i <- state$S_red_nodes[k]

        len_s_red <- length(state$S_red_nodes)
        state$S_red_nodes[k] <- state$S_red_nodes[len_s_red]
        state$S_red_nodes <- state$S_red_nodes[-len_s_red]

        state$I_red <- state$I_red + 1
        state$S_red  <- state$total_red  - state$I_red  - state$R_red
        state$I_red_nodes <- c(state$I_red_nodes, i)

      } else {

        # 1.2.  remove from blue susceptible set (swap-and-pop O(1))
        if (length(state$S_blue_nodes) == 0)
          return(state)

        k <- sample.int(length(state$S_blue_nodes), 1)
        i <- state$S_blue_nodes[k]

        len_s_blue <- length(state$S_blue_nodes)
        state$S_blue_nodes[k] <- state$S_blue_nodes[len_s_blue]
        state$S_blue_nodes <- state$S_blue_nodes[-len_s_blue]

        state$I_blue <- state$I_blue + 1
        state$S_blue <- state$total_blue - state$I_blue - state$R_blue
        state$I_blue_nodes <- c(state$I_blue_nodes, i)

      }

        # Infection history
        state$inf_time <- c(state$inf_time, t)
        state$inf_camp <- c(
          state$inf_camp,
          ifelse(state$opinions[i] < 0, "red", "blue")
      )

        state$S_count <- state$S_count - 1
        state$I_count <- state$I_count + 1
        state$I_nodes <- c(state$I_nodes, i)

       # Remove from S_nodes (swap-and-pop O(1))
        idx_s <- match(i, state$S_nodes)
        if (!is.na(idx_s)) {
          len_s <- length(state$S_nodes)
          state$S_nodes[idx_s] <- state$S_nodes[len_s]
          state$S_nodes <- state$S_nodes[-len_s]
        }


      state$epi[i] <- state$I

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

      state$I_count <- state$I_count - 1
      state$R_count <- state$R_count + 1

      # remove from I_nodes (swap-and-pop O(1))
      len_i_nodes <- length(state$I_nodes)
      state$I_nodes[k] <- state$I_nodes[len_i_nodes]
      state$I_nodes <- state$I_nodes[-len_i_nodes]

      # add to R_nodes
      state$R_nodes <- c(state$R_nodes, i)

      # update epidemic state
      state$epi[i] <- state$R

      # update camps
      if (state$opinions[i] < 0) {

        state$I_red <- state$I_red - 1
        state$R_red <- state$R_red + 1

        idx_ired <- match(i, state$I_red_nodes)
        if (!is.na(idx_ired)) {
          len_ired <- length(state$I_red_nodes)
          state$I_red_nodes[idx_ired] <- state$I_red_nodes[len_ired]
          state$I_red_nodes <- state$I_red_nodes[-len_ired]
        }

        state$R_red_nodes <-
          c(state$R_red_nodes, i)

      } else {

        state$I_blue <- state$I_blue - 1
        state$R_blue <- state$R_blue + 1

        idx_iblue <- match(i, state$I_blue_nodes)
        if (!is.na(idx_iblue)) {
          len_iblue <- length(state$I_blue_nodes)
          state$I_blue_nodes[idx_iblue] <- state$I_blue_nodes[len_iblue]
          state$I_blue_nodes <- state$I_blue_nodes[-len_iblue]
        }

        state$R_blue_nodes <-
          c(state$R_blue_nodes, i)

      }


    }


  return(state)
}

