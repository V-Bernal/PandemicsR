#' get_palette
#'
#' @docType package
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name get_palette
#'
#' @import Matrix
#' @import igraph
#'
#' @references \url{}
#' @seealso \code{\link{brocolors}}
#' @keywords hplot
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

  colorRampPalette(c("darkred", "lightcoral", "lightblue", "darkblue"))(num_opinions)

  }
