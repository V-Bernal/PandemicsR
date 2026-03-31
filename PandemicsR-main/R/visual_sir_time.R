#' Plot SIR Dynamics Over Time
#'
#' @author OpenAI Codex
#'
#' @name visual_sir_time
#'
#' @param sir_mat Matrix with columns \code{S}, \code{I}, and \code{R}.
#' @param time_points Optional vector of recorded simulation times.
#'
#' @return NULL
#' @export
visual_sir_time <- function(sir_mat, time_points = NULL) {
  sir_mat <- as.matrix(sir_mat)

  if (is.null(time_points) || length(time_points) != nrow(sir_mat)) {
    time_points <- seq_len(nrow(sir_mat)) - 1
    xlab <- "Recorded step"
  } else {
    time_points <- as.numeric(time_points)
    xlab <- "Simulation time"
  }

  xlim <- range(time_points)
  if (xlim[[1]] == xlim[[2]]) {
    xlim <- xlim + c(-0.5, 0.5)
  }

  ylim <- range(sir_mat)
  if (!all(is.finite(ylim)) || ylim[[1]] == ylim[[2]]) {
    ylim <- c(0, max(1, ylim[[2]]))
  }

  plot.new()
  plot.window(xlim = xlim, ylim = ylim)
  axis(1)
  axis(2)
  box()
  title(xlab = xlab, ylab = "Individuals")

  line_cols <- c(S = "#4d4d4d", I = "#b22222", R = "#1f77b4")
  for (label in colnames(sir_mat)) {
    lines(time_points, sir_mat[, label], lwd = 2, col = line_cols[[label]])
  }

  legend("right", legend = colnames(sir_mat), col = line_cols[colnames(sir_mat)], lwd = 2, bty = "n")
}
