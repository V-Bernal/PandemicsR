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
#' @param time_points optional vector of recorded simulation times
#'
#' @return NULL
#'
#' @examples
#' #---------------------------------
#' # visual_step_time
#' # visual_step_time(frac_mat, num_opinions)
#' #---------------------------------
#' @export
visual_step_time <- function(frac_mat, num_opinions, time_points = NULL){

  levels_vec <- get_levels_vec(num_opinions)
  my_palette <- get_palette(num_opinions)

  if (is.null(time_points) || length(time_points) != nrow(frac_mat)) {
    time_points <- seq_len(nrow(frac_mat)) - 1
    xlab <- "Recorded step"
  } else {
    time_points <- as.numeric(time_points)
    xlab <- "Simulation time"
  }

  xlim <- range(time_points)
  if (xlim[1] == xlim[2]) {
    xlim <- xlim + c(-0.5, 0.5)
  }

  plot.new()
  plot.window(xlim = xlim, ylim = c(0, 1))
  axis(1)
  axis(2)
  box()
  title(xlab = xlab, ylab = "Opinion fraction")

  for (k in 1:num_opinions) {
    lines(x = time_points, y = frac_mat[,k ],
          lwd = 2, col = my_palette[k])
  }

  # Update line plot for fraction of +1
  #Sys.sleep(0.05) # pause so you can see the update
}
