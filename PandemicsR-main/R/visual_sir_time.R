#' Plot SIR Dynamics Over Time
#'
#' @author OpenAI Codex
#'
#' @name visual_sir_time
#'
#' @param sir_mat Matrix with columns \code{S}, \code{I}, and \code{R}.
#' @param time_points Optional vector of recorded simulation times.
#' @param ylab Y-axis label.
#' @param main Optional plot title.
#'
#' @return NULL
#' @export
visual_sir_time <- function(sir_mat, time_points = NULL, ylab = "Individuals", main = NULL) {
  sir_mat <- as.matrix(sir_mat)
  if (is.null(colnames(sir_mat))) {
    colnames(sir_mat) <- c("S", "I", "R")[seq_len(ncol(sir_mat))]
  }

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

  line_cols <- c(S = "#4d4d4d", I = "#b22222", R = "#1f77b4")
  plot_cols <- unname(line_cols[colnames(sir_mat)])
  missing_cols <- is.na(plot_cols)
  if (any(missing_cols)) {
    plot_cols[missing_cols] <- grDevices::rainbow(sum(missing_cols))
  }
  names(plot_cols) <- colnames(sir_mat)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(5, 4.5, 5.5, 1.5))

  plot.new()
  plot.window(xlim = xlim, ylim = ylim, bty = "L")
  axis(1)
  axis(2)
  box()
  title(main = main, xlab = xlab, ylab = ylab)

  for (label in colnames(sir_mat)) {
    lines(time_points, sir_mat[, label], type = "s", lwd = 2, col = plot_cols[[label]])
  }

  legend(
    "top",
    legend = colnames(sir_mat),
    col = plot_cols[colnames(sir_mat)],
    lwd = 2,
    lty = 1,
    horiz = TRUE,
    bty = "n",
    inset = c(0, -0.18),
    xpd = NA
  )
}
