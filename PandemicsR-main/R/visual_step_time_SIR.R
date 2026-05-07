#' visual_step_time_SIR
#'
#' @author Victor Bernal
#'
#' @name visual_step_time_SIR
#'
#' @import graphics
#'
#' @param x SIR data frame
#'
#' @return NULL
#'
#' @examples
#' #---------------------------------
#' # visual_step_time_SIR(x)
#' #---------------------------------
#' @export
visual_step_time_SIR <- function(x){

  # increase top margin
  op <- par(mar = c(5, 4, 6, 2))  # bottom, left, top, right

  # actual simulation time
  time <- x[,1]

  plot.new()
  #plot.window(xlim = c(1, nrow(x)), ylim = c(0, 1), xaxt = "n", bty='L')  # suppress x-axis)
  plot.window(xlim = c(min(time), max(time)), ylim = c(0, 1), xaxt = "n", bty='L')  # suppress x-axis)

  # evenly spaced tick marks
  ticks <- pretty(time, n = 6)

  axis(1,
       at = ticks,
       labels = round(ticks, 2))

  #axis(1, at = c(1: nrow(x)), labels = (1:nrow(x)))
  axis(2)
  box()
  title(xlab = "Time", ylab = "Fraction SIR")

  cols <-  c("grey30", "yellow2", "seagreen")

  for (k in 2:ncol(x)) {
    lines(x = time, y = x[, k], type="s",
          lwd = 2, col = cols[k-1])
  }

  # legend ABOVE plot
  legend("top",
         legend = c("S","I","R"),
         col = cols,
         lwd = 2,
         lty = 1,
         horiz = TRUE,
         bty = "n",
         inset = c(0, -0.2),
         xpd = NA)

  # reset par
  par(op)

}
