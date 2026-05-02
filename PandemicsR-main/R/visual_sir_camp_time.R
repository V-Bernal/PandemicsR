#' Plot SIR Dynamics By Camp
#'
#' @author OpenAI Codex
#'
#' @name visual_sir_camp_time
#'
#' @param camp_sir_history List of camp-by-compartment SIR count matrices.
#' @param time_points Optional vector of recorded simulation times.
#'
#' @return NULL
#' @export
visual_sir_camp_time <- function(camp_sir_history, time_points = NULL) {
  if (is.list(camp_sir_history)) {
    if (!length(camp_sir_history)) {
      stop("camp_sir_history must contain at least one snapshot.")
    }
    camp_sir_array <- simplify2array(camp_sir_history)
  } else {
    camp_sir_array <- camp_sir_history
  }

  if (length(dim(camp_sir_array)) != 3L) {
    stop("camp_sir_history must be a list of matrices or a 3D array.")
  }

  camp_labels <- dimnames(camp_sir_array)[[1]]
  compartment_labels <- dimnames(camp_sir_array)[[2]]
  n_snapshots <- dim(camp_sir_array)[[3]]

  if (is.null(camp_labels)) {
    camp_labels <- paste("Camp", seq_len(dim(camp_sir_array)[[1]]))
  }
  if (is.null(compartment_labels)) {
    compartment_labels <- c("S", "I", "R")[seq_len(dim(camp_sir_array)[[2]])]
  }

  if (is.null(time_points) || length(time_points) != n_snapshots) {
    time_points <- seq_len(n_snapshots) - 1
    xlab <- "Recorded step"
  } else {
    time_points <- as.numeric(time_points)
    xlab <- "Simulation time"
  }

  xlim <- range(time_points)
  if (xlim[[1]] == xlim[[2]]) {
    xlim <- xlim + c(-0.5, 0.5)
  }

  line_cols <- c(S = "#4d4d4d", I = "#b22222", R = "#1f77b4")
  plot_cols <- unname(line_cols[compartment_labels])
  missing_cols <- is.na(plot_cols)
  if (any(missing_cols)) {
    plot_cols[missing_cols] <- grDevices::rainbow(sum(missing_cols))
  }
  names(plot_cols) <- compartment_labels

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(length(camp_labels), 1), mar = c(4.5, 4.5, 4.8, 1.5))

  for (camp_idx in seq_along(camp_labels)) {
    sir_mat <- t(camp_sir_array[camp_idx, , , drop = TRUE])
    colnames(sir_mat) <- compartment_labels
    ylim <- range(sir_mat)
    if (!all(is.finite(ylim)) || ylim[[1]] == ylim[[2]]) {
      ylim <- c(0, max(1, ylim[[2]]))
    }

    plot.new()
    plot.window(xlim = xlim, ylim = ylim, bty = "L")
    axis(1)
    axis(2)
    box()
    title(main = camp_labels[[camp_idx]], xlab = xlab, ylab = "Individuals")

    for (label in compartment_labels) {
      lines(time_points, sir_mat[, label], type = "s", lwd = 2, col = plot_cols[[label]])
    }

    legend(
      "top",
      legend = compartment_labels,
      col = plot_cols[compartment_labels],
      lwd = 2,
      lty = 1,
      horiz = TRUE,
      bty = "n",
      inset = c(0, -0.18),
      xpd = NA
    )
  }
}
