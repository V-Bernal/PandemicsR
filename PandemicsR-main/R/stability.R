#' Stability monitoring
#'
#' @export
# ------------------------------------------------------------
# Stability monitoring
# ------------------------------------------------------------

create_stability_monitor <- function(
    window = 100,
    threshold = 0.01,
    persistence = 5
) {

  list(
    window = window,
    threshold = threshold,
    persistence = persistence,

    history = list(),
    stable_windows = 0,
    is_stable = FALSE
  )
}


# ------------------------------------------------------------
# Update stability monitor
# ------------------------------------------------------------
update_stability <- function(
    monitor,
    state,
    time
) {

  metrics <- calculate_stability_metrics(state)

  monitor$history[[length(monitor$history) + 1]] <- list(
    time = time,
    metrics = metrics
  )


  if(time > monitor$persistence * monitor$window){
    monitor$history[1] <- NULL
  }

  stability <- calculate_stability(
    monitor$history,
    window = monitor$window,
    threshold = monitor$threshold
  )

  if (stability) {
    monitor$stable_windows <- monitor$stable_windows + 1
  } else {
    monitor$stable_windows <- 0
  }

  if (monitor$stable_windows >= monitor$persistence) {
    monitor$is_stable <- TRUE
  }

  return(monitor)
}
# ------------------------------------------------------------
# Generic metric extraction
# ------------------------------------------------------------

calculate_stability_metrics <- function(state) {

  # Placeholder.
  #
  # Later we can plug in:
  #   - mean opinion
  #   - opinion variance
  #   - number of groups
  #   - largest group
  #   - group-size distribution
  #   - event rates
  #   - etc.

  list(
    mean_opinion = mean(state$opinions),
    n_groups = length(state$members)
  )
}


# ------------------------------------------------------------
# Generic stability calculation
# ------------------------------------------------------------

calculate_stability <- function(
    history,
    window,
    threshold
) {
  # Placeholder.
  #
  # Later this function decides what
  # "stable" actually means.


    if (length(history) < 2) {
      return(FALSE)
    }

    times <- sapply(history, function(x) x$time)

    current <- history[[length(history)]]

    previous_index <- which(
      times <= current$time - window
    )

    if (length(previous_index) == 0) {
      return(FALSE)
    }

    previous <- history[[max(previous_index)]]

    opinion_change <- abs(
      current$metrics$mean_opinion -
        previous$metrics$mean_opinion
    )

    groups_unchanged <-
      current$metrics$n_groups ==
      previous$metrics$n_groups

    opinion_change < threshold &&
      groups_unchanged
  }
