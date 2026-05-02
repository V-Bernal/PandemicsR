#' Plot Final Opinion Shares
#'
#' @author OpenAI Codex
#'
#' @name visual_opinion_shares
#'
#' @param frac_mat Matrix of opinion fractions over time.
#' @param num_opinions Number of opinion states.
#'
#' @return NULL
#' @export
visual_opinion_shares <- function(frac_mat, num_opinions) {
  frac_mat <- as.matrix(frac_mat)
  final_share <- frac_mat[nrow(frac_mat), ]
  if (is.null(names(final_share))) {
    names(final_share) <- get_state_labels(num_opinions)
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(6, 4.5, 4, 1))

  bar_cols <- get_palette(num_opinions)
  barplot(
    height = final_share,
    names.arg = names(final_share),
    col = bar_cols[seq_along(final_share)],
    ylim = c(0, 1),
    ylab = "Share",
    main = "Final opinion shares",
    las = 2
  )
}
