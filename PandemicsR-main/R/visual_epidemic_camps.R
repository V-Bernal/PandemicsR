#' Plot Final Epidemic Outcomes By Camp
#'
#' @author OpenAI Codex
#'
#' @name visual_epidemic_camps
#'
#' @param camp_attack_rate Named numeric vector of final attack rates.
#'
#' @return NULL
#' @export
visual_epidemic_camps <- function(camp_attack_rate) {
  camp_attack_rate <- as.numeric(camp_attack_rate)
  label_names <- names(camp_attack_rate)
  if (is.null(label_names) || any(!nzchar(label_names))) {
    label_names <- get_camp_labels()
  }
  names(camp_attack_rate) <- label_names

  bar_cols <- c("#b22222", "#1f77b4")
  bp <- barplot(
    height = camp_attack_rate,
    names.arg = names(camp_attack_rate),
    col = bar_cols[seq_along(camp_attack_rate)],
    ylim = c(0, 1),
    ylab = "Attack rate",
    main = "Final epidemic outcomes by camp"
  )
  text(bp, camp_attack_rate, labels = sprintf("%.2f", camp_attack_rate), pos = 3, cex = 0.9)
}
