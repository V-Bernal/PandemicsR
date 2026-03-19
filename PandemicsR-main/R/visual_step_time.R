#' visual_step_time
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name visual_step_time
#'
#' @import graphics
#'
#' @param frac_mat fraction of opinions
#' @param num_opinions number of opinions
#' @param n number of individuals
#'
#' @return NULL
#'
#' @examples
#' #---------------------------------
#' # visual_step_time
#' # visual_step_time(frac_mat, num_opinions)
#' #---------------------------------
#' @export
visual_step_time <- function(frac_mat, num_opinions, n){

  levels_vec <- get_levels_vec(num_opinions)
  my_palette <- get_palette(num_opinions)

  # Note: frac was recorded every n events. ticks are rescaled by mult by n
  plot.new()
  plot.window(xlim = c(1, nrow(frac_mat)), ylim = c(0, 1), xaxt = "n")  # suppress x-axis)
  axis(1, at = c(1: nrow(frac_mat)), labels = n*(1:nrow(frac_mat)))
  #plot.window(xlim = c(1, nrow(frac_mat)), ylim = c(0, 1))
  #axis(1)
  axis(2)
  box()
  title(xlab = "Time", ylab = "Voter Fraction")

  for (k in 1:num_opinions) {
    lines(x = 1:nrow(frac_mat), y = frac_mat[,k ],
          lwd = 2, col = my_palette[k])
  }

  # Update line plot for fraction of +1
  #Sys.sleep(0.05) # pause so you can see the update
}
