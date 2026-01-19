#' visual_degree_time
#'
#' @docType package
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name visual_degree_time
#'
#' @param avg_deg avg_deg
#'
#' @return NULL
#'
#' @examples
#' #---------------------------------
#' # visual_degree_time
#' # visual_degree_time(avg_deg)
#' #---------------------------------
#' @export
visual_degree_time <- function(avg_deg){

  plot.new()
  plot.window(xlim = c(1, length(avg_deg)), ylim = c(0, 1))
  axis(1)
  axis(2)
  box()

  lines(x = 1:length(avg_deg), y = avg_deg, lwd = 2)

  # Update line plot for fraction of +1
  #Sys.sleep(0.05) # pause so you can see the update
}


