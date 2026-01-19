#' Single Schelling step
#'
#' @docType package
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name schelling_step
#'
#' @import Matrix
#' @import igraph
#'
#' @param RIG number of individuals vertices
#' @param opinions number of group vertices
#' @param threshold threshold
#'
#' @return rig
#'
#' @examples
#' #---------------------------------
#' # A single Schelling step
#' # schelling_step(RIG, opinions, threshold)
#' #---------------------------------
#' @export

schelling_step <- function(RIG, opinions, threshold){
  n <- length(opinions)
  for(i in 1:n){
    neighbors <- which(RIG[i,] == 1)
    if(length(neighbors) == 0) next
    s_i <- sum(opinions[i] == opinions[neighbors])/length(neighbors)
    if(s_i < threshold){
      # remove one dissimilar neighbor
      dissimilar <- neighbors[opinions[neighbors] != opinions[i]]
      if(length(dissimilar) > 0){
        remove_j <- sample(dissimilar,1)
        RIG[i, remove_j] <- 0
        RIG[remove_j, i] <- 0
      }
      # add edge to a similar node not already connected
      candidates <- which(opinions == opinions[i] & RIG[i,] == 0 & (1:n) != i)
      if(length(candidates) > 0){
        add_j <- sample(candidates,1)
        RIG[i, add_j] <- 1
        RIG[add_j, i] <- 1
      }
    }
  }
  return(RIG)
}
