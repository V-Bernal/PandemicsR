#' Initialize voters' opinions
#'
#' @docType package
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name initialize_opinions
#'
#' @import Matrix
#' @import igraph
#'
#' @references \url{}
#' @seealso \code{\link{brocolors}}
#' @keywords hplot
#'
#' @param n number of individuals vertices
#'
#' @return initialize_opinions \code{x}.
#'
#' @examples
#' #---------------------------------
#' # Initialize opinions
#' # initialize_opinions(n = 20)
#' #---------------------------------
#' @export

initialize_opinions <- function(n) {
  sample(c(-1, 1), n, replace = TRUE)
}
