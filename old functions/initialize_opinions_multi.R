#' Initialize voters' opinions multi
#'
#' @docType package
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name initialize_opinions_multi
#'
#' @import Matrix
#' @import igraph
#'
#' @references \url{}
#' @seealso \code{\link{brocolors}}
#' @keywords hplot
#'
#' @param n number of individuals vertices
#' @param num_opinions number of opinions
#'
#' @return initialize_opinions \code{x}.
#'
#' @examples
#' #---------------------------------
#' # Initialize opinions
#' # initialize_opinions_multi(n = 20,num_Opinions =2)
#' #---------------------------------
#' @export

initialize_opinions_multi <- function(n, num_opinions = NULL) {

  set_opinions <- get_levels_vec(num_opinions)

  opinions_sample <- sample(set_opinions , n, replace = TRUE)

  return(opinions_sample)

  }
