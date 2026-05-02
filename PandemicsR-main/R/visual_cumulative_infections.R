#' Plot Cumulative Infections By Camp
#'
#' @author OpenAI Codex
#'
#' @name visual_cumulative_infections
#'
#' @param infection_events Data frame with \code{time} and \code{camp} columns.
#'   An optional \code{count} column is used for aggregate simulations.
#' @param final_time Optional final simulation time for the x-axis.
#' @param camp_labels Labels for the camps to draw.
#'
#' @return NULL
#' @export
visual_cumulative_infections <- function(infection_events, final_time = NULL, camp_labels = get_camp_labels()) {
  if (is.null(infection_events)) {
    infection_events <- data.frame(time = numeric(0), camp = character(0))
  }

  infection_events <- as.data.frame(infection_events)
  if (!all(c("time", "camp") %in% names(infection_events))) {
    stop("infection_events must contain time and camp columns.")
  }

  event_times <- as.numeric(infection_events$time)
  event_camps <- as.character(infection_events$camp)
  event_counts <- if ("count" %in% names(infection_events)) {
    as.numeric(infection_events$count)
  } else {
    rep(1, length(event_times))
  }
  finite_idx <- is.finite(event_times) & is.finite(event_counts)
  finite_times <- event_times[finite_idx]

  x_max <- max(c(1, final_time, finite_times), na.rm = TRUE)
  y_max <- max(1, vapply(camp_labels, function(camp_label) {
    sum(event_counts[finite_idx & event_camps == camp_label])
  }, numeric(1)))
  xlim <- c(0, x_max)
  ylim <- c(0, y_max)

  line_cols <- c("#b22222", "#1f77b4", grDevices::rainbow(max(0, length(camp_labels) - 2L)))
  line_cols <- line_cols[seq_along(camp_labels)]
  names(line_cols) <- camp_labels

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(5, 4.5, 5.5, 1.5))

  plot.new()
  plot.window(xlim = xlim, ylim = ylim, bty = "L")
  axis(1)
  axis(2)
  box()
  title(xlab = "Simulation time", ylab = "Cumulative infections")

  for (camp_label in camp_labels) {
    camp_idx <- which(event_camps == camp_label & is.finite(event_times) & is.finite(event_counts))
    camp_idx <- camp_idx[order(event_times[camp_idx])]
    camp_times <- event_times[camp_idx]
    camp_counts <- event_counts[camp_idx]
    y_values <- cumsum(camp_counts)
    if (!length(camp_times)) {
      lines(c(0, x_max), c(0, 0), type = "s", lwd = 2, col = line_cols[[camp_label]])
      next
    }

    lines(
      c(0, camp_times, x_max),
      c(0, y_values, tail(y_values, 1)),
      type = "s",
      lwd = 2,
      col = line_cols[[camp_label]]
    )
  }

  legend(
    "top",
    legend = camp_labels,
    col = line_cols[camp_labels],
    lwd = 2,
    lty = 1,
    horiz = TRUE,
    bty = "n",
    inset = c(0, -0.18),
    xpd = NA
  )
}
