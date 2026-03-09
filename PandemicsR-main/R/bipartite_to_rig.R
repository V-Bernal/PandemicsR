#' Bipartite to RIG
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name bipartite_to_rig
#'
#' @import Matrix
#' @import igraph
#'
#'
#' @param bipartite_matrix bipartite_matrix
#'
#' @return Random Intersection Graph \code{rig}.
#'
#' @examples
#' #---------------------------------
#' # bipartite_to_rig(bipartite_matrix)
#' #---------------------------------
#' @export
bipartite_to_rig <- function(bipartite_matrix) {
  rig <- bipartite_matrix %*% t(bipartite_matrix)
  rig@x[rig@x > 0] <- 1
  diag(rig) <- 0 # avoid self loops
  return(rig)
}
