#' Single Schelling step_ on the bipartite graph with vertex at random
#'
#' @docType package
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name schelling_step_bipartite_random
#'
#' @import Matrix
#' @import igraph
#'
#' @param B bipartite grpah
#' @param opinions number of group vertices
#' @param threshold threshold
#' @param vertex vertex index
#' @param RIG RIG
#' @param const cont
#'
#' @return bipartite_matrix, RIG
#'
#' @examples
#' #---------------------------------
#' # A single Schelling step
#' # schelling_step_bipartite_random(bipartite_matrix, opinions, vertex)
#' #---------------------------------
#' @export
schelling_step_bipartite_random <- function(bipartite_matrix, RIG, opinions, vertex, const = 1) {
  #-------------------
  # Description: This function runs a Schelling edges' update over the aux bipartite graph.
  # The bipartite graph has n individuals and m groups
  # For vertex i (indivdual), loop over all the m groups.
  # Change edges (addition and deletion) with prob beta
  # Prob beta is increasing in the percentage of other group members opinions
  # Watch out*: Why - alpha? Denominator can be zero?
  #-------------------

  # node i
  i <- vertex

  n <- nrow(bipartite_matrix)
  m <- ncol(bipartite_matrix)

  # beta
  #const <- 1
  p_addition <- const / m
  alpha <- opinions[i]

  # For each neighbor m, try addition and deletion
    for (j in 1:m) {

      k_pos <- sum(bipartite_matrix[, j] * (opinions == 1)) # edges with opinion 1 including i itself
      k_neg <- sum(bipartite_matrix[, j] * (opinions == -1)) # edges with opinion -1 including i itself

      new_value <- bipartite_matrix[i, j] # in case of no update

      # Probability of edge i - j creation
      if (bipartite_matrix[i, j] == 0) {
        if (runif(1) < p_addition) {
          #print(p_addition)
          new_value <- 1
        }
      }

      # Probability of edge i - j deletion
      if (bipartite_matrix[i, j] == 1) {
        p_deletion <- -1 * alpha * (k_pos - k_neg) / (k_pos + k_neg) # Watch out*: added +1 deno zero?
        # p_deletion <- -1 * alpha * (k_pos - k_neg) / (k_pos + k_neg) > 0.5
        #p_deletion <- 1- p_addition
        #p_deletion <- abs(alpha * (k_pos - k_neg) / (k_pos + k_neg) ) # Watch out*: added +1 deno zero?
        #print(p_deletion)
        if (runif(1) < p_deletion) {
          new_value <- 0
          }
      }

      # update bipartite
      bipartite_matrix[i, j] <- new_value

    }

  RIG <- bipartite_to_rig(bipartite_matrix)

  #----------------------------------
  # More efficient to update blocks, not the entire RIG.
  #----------------------------------
  # new_neighbors <- setdiff(as.numeric(which((bipartite_matrix[i, ] %*% t(bipartite_matrix)) > 0)), i)
  # neighbors_of_i <- unique(c(which(RIG[i, ] == 1), new_neighbors))
  # affected <- sort(unique(c(i, neighbors_of_i)))
  #
  # # Update only if there was a change (isolated node i has no changes)
  # if (length(affected) > 1) {
  #   sub_RIG <- Matrix((bipartite_matrix[affected, ] %*% t(bipartite_matrix[affected, ])) > 0, sparse = TRUE)
  #   diag(sub_RIG) <- 0
  #   RIG[affected, affected] <- Matrix(sub_RIG, sparse = TRUE)
  # }
  #----------------------------------
  return(list(bipartite_matrix, RIG))
}
