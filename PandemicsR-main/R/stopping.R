#' Stopping
#'
#' Determine whether the simulation should stop.
#'
#' @param params Simulation parameters.
#' @param state Current simulation state.
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

