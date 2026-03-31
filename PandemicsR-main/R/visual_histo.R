#' visual_histo
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name visual_histo
#'
#' @import graphics
#'
#' @param opinion_history opinion_history
#' @param num_opinions num_opinions
#'
#' @return NULL
#'
#' @examples
#' #---------------------------------
#' # visual_histo
#' # visual_histo(opinion_history, num_opinions)
#' #---------------------------------
#' @export
visual_histo <- function(opinion_history, num_opinions){

  levels_vec <- get_levels_vec(num_opinions)
  state_labels <- get_state_labels(num_opinions)
  my_palette <- get_palette(num_opinions)

  tb0 <- table(factor(opinion_history[,1], levels = levels_vec))
  tb <- table(factor(opinion_history[, ncol(opinion_history)], levels = levels_vec))
  m <- max(c(tb0, tb))

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(mfrow = c(1, 2))
  barplot(tb0, names.arg = state_labels, col = my_palette, ylim = c(0, m), las = 2, main = "Initial")
  barplot(tb, names.arg = state_labels, col = my_palette, ylim = c(0, m), las = 2, main = "Final")
}
