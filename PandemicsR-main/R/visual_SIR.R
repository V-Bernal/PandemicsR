#' visual_SIR
#'
#' @author Victor Bernal
#'
#' @name visual_SIR
#'
#' @import graphics
#'
#' @param x SIR data frame
#'
#' @return NULL
#'
#' @examples
#' #---------------------------------
#' # visual_SIR(x)
#' #---------------------------------
#' @export
visual_SIR <- function(inf_time, inf_camp){

  # time of infection
  red_times  <- inf_time[inf_camp == "red"]
  blue_times <- inf_time[inf_camp == "blue"]

  plot.new()
  #par(mar = c(5, 4, 4, 2))  # bottom, left, top, right
  plot.window(xlim = c(0, max(inf_time)), ylim = c(0, max(seq_along(inf_time))),
              xaxt = "n", bty='L')  # suppress x-axis)
  axis(1, at = c(1: max(inf_time)), labels = (1:max(inf_time)) )
  axis(2)
  box()
  title(xlab = "Time", ylab = "Cumulative infections")


  lines(x = sort(red_times), y = seq_along(red_times), lwd = 2, type="b", col = 'red')
  lines(x = sort(blue_times), y = seq_along(blue_times), lwd = 2, type="b", col = 'blue')#,
       #xlab="Time", ylab="Cumulative infections")
  #legend("topleft", legend=c("Red","Blue"), lwd=2,col = c("Blue", "Red"))
  }
