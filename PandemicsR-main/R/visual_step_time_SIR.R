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

  plot.new()
  #par(mar = c(5, 4, 4, 2))  # bottom, left, top, right
  plot.window(xlim = c(1, nrow(x)), ylim = c(0, 1), xaxt = "n", bty='L')  # suppress x-axis)
  axis(1, at = c(1: nrow(x)), labels = (1:nrow(x)))
  axis(2)
  box()
  title(xlab = "Time", ylab = "Fraction SIR")

  cols <-  c("grey30", "yellow2", "seagreen")

  for (k in 2:ncol(x)) {
    lines(x = 1:nrow(x), y = x[, k], type="s",
          lwd = 2, col = cols[k-1])
  }

  # legend(x = nrow(x)/2 , y = 1.2, col = cols, pch = 20, horiz = TRUE,
  #        legend=c("S","I","R"), y.intersp = 0.5, bty = "n",
  #        lty=1, lwd=2)

  # legend("top",
  #        legend = c("S","I","R"),lty=1,
  #        lwd = 2, col = cols, y.intersp = 0.35,
  #        inset = c(0, -0.15),  # move upward outside the box
  #        xpd = NA,   # allow drawing outside plot
  #        horiz = TRUE)         # optional: horizontal layout


}
