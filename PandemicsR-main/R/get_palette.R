#' get_palette
#'
#' @author Victor Bernal
#'
#' @name get_palette
#'
#' @param num_opinions number of opinions
#'
#' @return palette \code{x}.
#'
#' @examples
#' #---------------------------------
#' # Initialize color palette
#' # get_palette (num_opinions)
#' #---------------------------------
#' @export
get_palette <- function(num_opinions){

  grDevices::colorRampPalette(c("darkred", "lightcoral", "lightblue", "darkblue"))(num_opinions)

}
