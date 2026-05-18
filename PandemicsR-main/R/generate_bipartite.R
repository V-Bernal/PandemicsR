#' Generate a Bipartite Matrix
#'
#' Simulates a bipartite network between two sets of vertices (e.g.,
#' individuals and groups) using a Poisson edge-generation mechanism.
#' The expected number of edges between each pair of nodes is determined
#' by their associated weights.
#'
#' Each entry in the resulting matrix represents binary membership between
#' an individual and a group. Poisson draws greater than zero are collapsed
#' to one.
#'
#' @author Victor Bernal
#'
#' @name generate_bipartite
#'
#' @importFrom stats rpois
#'
#' @param n Integer. Number of vertices in the first set (e.g., individuals).
#' @param m Integer. Number of vertices in the second set (e.g., groups).
#' @param individual_weights Numeric vector of length \code{n}. Weights
#' associated with individual vertices, controlling their expected degree.
#' @param group_weights Numeric vector of length \code{m}. Weights associated
#' with group vertices, controlling their attractiveness.
#'
#' @return A sparse matrix of class \code{dgCMatrix} with dimension
#' \code{n x m}, representing the bipartite adjacency matrix. Entries are
#' binary, with one indicating membership.
#'
#' @details
#' The expected number of edges between individual \eqn{i} and group \eqn{j}
#' is given by:
#' \deqn{\lambda_{ij} = \frac{w_i \cdot g_j}{\sum_k g_k}}
#' where \eqn{w_i} and \eqn{g_j} are the individual and group weights,
#' respectively. Each entry is then sampled and binarized as:
#' \deqn{B_{ij} \sim \mathrm{Poisson}(\lambda_{ij})}
#'
#' This formulation produces a binary random bipartite graph consistent with
#' weighted random intersection membership models.
#' @examples
#' #---------------------------------
#' # Generate a bipartite
#' # Set parameters
#' n <- 5
#' m <- 3
#' individual_weights <- runif(n, 0.5, 1.5)
#' group_weights <- runif(m, 0.5, 1.5)
#'
#' # Generate bipartite matrix
#' B <- generate_bipartite(n, m, individual_weights, group_weights)
# Inspect result
#' print(B)
#' #---------------------------------
#' @export
generate_bipartite <- function(n, m, individual_weights, group_weights) {
  total_group_weight <- sum(group_weights)
  if (!is.finite(total_group_weight) || total_group_weight <= 0) {
    return(Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      x = integer(0),
      dims = c(n, m)
    ))
  }

  # Expected Poisson edges
  expected_edges <- tcrossprod(individual_weights, group_weights) / total_group_weight
  # Generate Poisson edges and collapse multi-edges to binary membership.
  edge_counts <- rpois(length(expected_edges), lambda = c(expected_edges))
  edge_counts[edge_counts > 0L] <- 1L

  bipartite_matrix <- Matrix::Matrix(
    matrix(
      edge_counts,
      nrow = n,
      ncol = m
    ),
    sparse = TRUE
  )
  return(bipartite_matrix)
}

