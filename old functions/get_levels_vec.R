#' get_levels_vec
#'
#' @docType package
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name get_levels_vec
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
#' @return levels
#'
#' @examples
#' #---------------------------------
#' # get_levels_vec (num_opinions)
#' #---------------------------------
#' @export
get_levels_vec <- function(num_opinions){

  lv <- seq(-num_opinions/2, num_opinions/2, by = 1)

  lv <- setdiff(lv, 0)

  return(lv)
}
