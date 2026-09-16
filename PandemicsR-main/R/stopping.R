#' Stopping
#'
#' Determine whether the simulation should stop.
#'
#' @param params Simulation parameters.
#' @param opinions Agent opinions.
#' @param B0 Initial network adjacency matrix.
#' @param state Current simulation state.
#' @param trackers Simulation trackers.
#' @param network_state Current network state.
#' @param stop_reason Current stopping reason.
#' @param t Current simulation time.
#' @param comp_time Computed event time.
#' @param event_counter Event counter.
#'
#' @return Updated stopping state.
#' @export
check_stopping <- function(state, params) {

  if (isTRUE(state$epidemic_started) && isTRUE(params$runEpidemic)) {

    if (state$I_count == 0) {
      return(list(
        stop = TRUE,
        reason = "No individuals are infected (epidemic ended)"
      ))
    }

    if (state$R_count == params$n) {
      return(list(
        stop = TRUE,
        reason = "All individuals recovered (epidemic ended)"
      ))
    }
  }

  if (params$num_opinions == 4) {

    if (all(state$opinions %in% c(-2, 2))) {
      return(list(
        stop = TRUE,
        reason = "Full polarization (all extreme opinions)"
      ))
    }
  }

  list(
    stop = FALSE,
    reason = NULL
  )
}

