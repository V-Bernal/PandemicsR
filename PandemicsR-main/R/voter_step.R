#' Single voter step
#'
#' @docType package
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name voter_step
#'
#' @import Matrix
#' @import igraph
#'
#' @param rig number of individuals vertices
#' @param opinions number of group vertices
#' @param kappa Poisson clock ringing rate (see \code{\link[graphics]{par}})
#' @param vertex vertex index
#'
#' @return opinions \code{opinions}.
#'
#' @examples
#' #---------------------------------
#' # A single voter step
#' # voter_step(rig, opinions, kappa)
#' #---------------------------------
#' @export
voter_step <- function(rig, opinions, kappa, vertex) {
  #-------------------
  # Description: This function runs a Voter opinion update.
  # For vertex i, choose a neighbor at random.
  # Change the opinion of vertex i to the opinion of neighbor j with prob kappa.
  # Watch out*: if(length(neighbors) > 1) because sample still samples 1:12.
  #-------------------

  # node i
  i <- vertex
  # i <- sample(1:length(opinions), 1)

  # node j
  if (runif(1) < kappa) {
    neighbors <- which(rig[i, ] != 0)

    if (length(neighbors) > 1) { # Watch out*
      j <- sample(x = neighbors, size = 1, replace = FALSE)
      opinions[i] <- opinions[j]
      #new_opinions <- opinions
      #new_opinions[i] <- opinions[j]
      #opinions <- new_opinions
    } else if(length(neighbors) == 1) { # Watch out*
      j <- neighbors
      opinions[i] <- opinions[j]
      #new_opinions <- opinions
      #new_opinions[i] <- opinions[j]
      #opinions <- new_opinions
    }

  }

  return(opinions)
}
