#' Plot Final Epidemic Outcomes By Camp
#'
#' @author OpenAI Codex
#'
#' @name visual_epidemic_camps
#'
#' @param final_camp_sir Matrix of final SIR counts by camp.
#'
#' @return NULL
#' @export
visual_epidemic_camps <- function(final_camp_sir) {
  final_camp_sir <- as.matrix(final_camp_sir)
  if (is.null(rownames(final_camp_sir))) {
    rownames(final_camp_sir) <- get_camp_labels()
  }

  camp_sizes <- rowSums(final_camp_sir)
  attack_rate <- ifelse(camp_sizes > 0, 1 - final_camp_sir[, "S"] / camp_sizes, NA_real_)
  plot_values <- ifelse(is.na(attack_rate), 0, attack_rate)
  text_labels <- ifelse(
    is.na(attack_rate),
    "No members",
    sprintf("%.2f", attack_rate)
  )
  label_names <- sprintf(
    "%s\nn=%s",
    rownames(final_camp_sir),
    format(round(camp_sizes), big.mark = ",", trim = TRUE)
  )
  bar_cols <- ifelse(camp_sizes > 0, c("#b22222", "#1f77b4")[seq_len(nrow(final_camp_sir))], "#d9d9d9")

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(6, 4.5, 4.5, 1))

  bp <- barplot(
    height = plot_values,
    names.arg = label_names,
    col = bar_cols,
    ylim = c(0, 1),
    ylab = "Attack rate",
    main = "Final epidemic outcomes by camp",
    las = 1
  )
  text(bp, pmax(plot_values, 0.02), labels = text_labels, pos = 3, cex = 0.9)
  mtext("Labels include final camp size. Empty camps are shown explicitly.", side = 1, line = 4.3, cex = 0.85, col = "gray35")
}
