#' heatmapPlot
#'
#' @docType package
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name heatmapPlot
#'
#'
#' @param opinion_history opinion_history
#' @param num_opinions num_opinions
#'
#' @return NULL
#'
#' @examples
#' #---------------------------------
#' # heatmapPlot
#' # heatmapPlot(opinion_history, num_opinions)
#' #---------------------------------
#' @export
heatmapPlot <- function(opinion_history, num_opinions){

  levels_vec <- get_levels_vec(num_opinions)
  my_palette <- get_palette(num_opinions)

  # Map opinions to numeric indices for colors
  mat_indices <- match(opinion_history, levels_vec)  # matrix of 1..num_opinions
  image(
    t(opinion_history),
    col = my_palette,
    axes = FALSE,
    xlab = "Time",
    ylab = "Individuals",
    main = "Opinion Evolution"
  )
  grid(ny = nrow(opinion_history), nx =  NA)

  #axis(1, at = seq(0, 1, length.out = ncol(mat)), labels = '')
  #axis(2, at = seq(0, 1, length.out = nrow(mat)), labels = nrow(mat):1)
}
