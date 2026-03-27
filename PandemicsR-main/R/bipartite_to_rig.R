#' Convert a Bipartite Matrix to a Random Intersection Graph
#'
#' Computes the one-mode projection of a bipartite network onto the set of
#' row nodes, producing a Random Intersection Graph (RIG). Two nodes are
#' connected in the resulting graph if they share at least one common
#' neighbor in the bipartite structure.
#'
#' The function performs a matrix multiplication of the bipartite matrix
#' with its transpose and then binarizes the result. Diagonal elements are
#' set to zero to avoid self-loops.
#'
#' @author Victor Bernal
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
