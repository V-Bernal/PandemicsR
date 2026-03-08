#' Generate bipartite matrix
#'
#' @docType package
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name generate_bipartite
#'
#' @import Matrix
#' @import igraph
#'
#' @references \url{}
#' @seealso \code{\link{brocolors}}
#' @keywords hplot
#'
#' @param n number of individuals vertices
#' @param m number of group vertices
#' @param individual_weights ind vertex weight
#' @param group_weights (see \code{\link[graphics]{par}})
#'
#' @return bipartite matrix \code{x}.
#'
#' @examples
#' #---------------------------------
#' # Generate a bipartite
#' # generate_bipartite(n = 20, m = 4, individual_weights, group_weights, lambda)
#' #---------------------------------
#' @export
#'
# generate_bipartite <- function(n, m, individual_weights, group_weights, lambda) {
#   expected_edges <- outer(individual_weights, group_weights,
#                           FUN = function(x, y) (x * y) / sum(group_weights))
#   bipartite_matrix <- Matrix::Matrix(rpois(n * m, lambda = lambda * expected_edges) > 0,
#                              nrow = n, ncol = m, sparse = TRUE)
#   diag(bipartite_matrix) <- 0
#   return(bipartite_matrix)
# }

generate_bipartite <- function(n, m, individual_weights, group_weights) {
  # Expected Poisson edges
  expected_edges <- outer(individual_weights, group_weights,
                          FUN = function(x, y) (x * y) / sum(group_weights))
  # Generate Poisson multi-edges
  bipartite_matrix <- Matrix::Matrix(rpois(n * m, lambda = expected_edges),
                                     nrow = n, ncol = m, sparse = TRUE)
  return(bipartite_matrix)
}






