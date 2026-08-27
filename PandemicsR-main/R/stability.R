#' Stability monitoring
#' @param monitor
#' @param state
#' @param time
#' @export
# ------------------------------------------------------------
# Stability monitoring
# ------------------------------------------------------------
create_stability_monitor <- function(
   window = stability_window,
    threshold = stability_threshold
) {

  list(
    window = window,
    threshold = threshold,

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

  # Add current observation
  monitor$history[[length(monitor$history) + 1]] <- list(
    time = time,
    metrics = metrics
  )

  # Remove observations that are no longer needed
  cutoff <- time - monitor$window

  monitor$history <- Filter(
    function(x) x$time >= cutoff,
    monitor$history
  )

  # Calculate stability over the window
  monitor$is_stable <- calculate_stability(
    monitor$history,
    window = monitor$window,
    threshold = monitor$threshold
  )

  monitor
}
# ------------------------------------------------------------
# Generic metric extraction
# ------------------------------------------------------------
calculate_stability_metrics <- function(state) {

  list(
    mean_opinion = mean(state$opinions),
    n_groups = length(state$members)
  )
}

#==============
calculate_stability <- function(history, window, threshold) {
  TRUE
  # if (length(history) < 2)
  #   return(FALSE)
  #
  # current <- history[[length(history)]]
  #
  # times <- sapply(history, function(x) x$time)
  #
  # previous_index <- which(
  #   times <= current$time - window
  # )
  #
  # if (length(previous_index) == 0)
  #   return(FALSE)
  #
  # previous <- history[[max(previous_index)]]
  #
  # opinion_change <- abs(
  #   current$metrics$mean_opinion -
  #     previous$metrics$mean_opinion
  # )
  #
  # groups_unchanged <-
  #   current$metrics$n_groups ==
  #   previous$metrics$n_groups
  #
  # opinion_change < threshold &&
  #   groups_unchanged
}
