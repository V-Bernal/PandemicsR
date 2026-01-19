#' Single Schelling step bipartite
#'
#' @docType package
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name schelling_step_bipartite
#'
#' @import Matrix
#' @import igraph
#'
#' @param B bipartite grpah
#' @param opinions number of group vertices
#' @param threshold threshold
#'
#' @return B
#'
#' @examples
#' #---------------------------------
#' # A single Schelling step
#' # schelling_step(B, opinions, threshold)
#' #---------------------------------
#' @export

schelling_step_bipartite <- function(B, opinions, threshold){
  n <- nrow(B)
  m <- ncol(B)

  for(i in 1:n){
    # Groups the individual belongs to
    groups_i <- which(B[i, ] == 1)

    if(length(groups_i) == 0) next  # skip isolated nodes

    # Compute satisfaction: fraction of neighbors in same opinion across groups
    same_count <- 0
    total_neighbors <- 0
    for(g in groups_i){
      members <- which(B[, g] == 1)        # all individuals in group g
      same_count <- same_count + sum(opinions[members] == opinions[i]) - 1  # exclude self
      total_neighbors <- total_neighbors + length(members) - 1
    }
    satisfaction_i <- if(total_neighbors>0) same_count / total_neighbors else 1

    # Rewire if below threshold
    if(satisfaction_i < threshold){

      # --- Remove edge from a group with many dissimilar neighbors ---
      dissimilar_groups <- groups_i[sapply(groups_i, function(g){
        members <- which(B[,g]==1)
        mean(opinions[members] != opinions[i]) > 0
      })]
      if(length(dissimilar_groups) > 0){
        remove_g <- sample(dissimilar_groups, 1)
        B[i, remove_g] <- 0
      }

      # --- Add edge to a group with more similar neighbors, not yet connected ---
      candidate_groups <- setdiff(1:m, which(B[i,]==1))
      similar_groups <- candidate_groups[sapply(candidate_groups, function(g){
        members <- which(B[,g]==1)
        if(length(members)==0) return(TRUE)  # empty group
        mean(opinions[members] == opinions[i]) > 0
      })]
      if(length(similar_groups) > 0){
        add_g <- sample(similar_groups, 1)
        B[i, add_g] <- 1
      }
    }
  }

  # Return updated bipartite matrix
  return(B)
}

