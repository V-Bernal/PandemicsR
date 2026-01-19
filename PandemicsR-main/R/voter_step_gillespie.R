#' Single voter_step_gillespie
#'
#' @docType package
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name voter_step_gillespie
#'
#' @import Matrix
#' @import igraph
#'
#' @param rig number of individuals vertices
#' @param opinions number of group vertices
#' @param kappa Poisson clock ringing rate (see \code{\link[graphics]{par}})
#'
#' @return opinions \code{opinions}.
#'
#' @examples
#' #---------------------------------
#' # voter_step_gillespie
#' # voter_step_gillespie(rig, opinions, kappa)
#' #---------------------------------
#' @export
voter_step_gillespie <- function(rig, opinions, kappa) {

  # Parameters
  N <- length(opinions)

  # 1. Compute rates for all nodes
  rates <- sapply(1:N, function(i) {
    if (sum(rig[i, ] != 0) > 0) kappa else 0
  })

  R_tot <- sum(rates)
  if (R_tot == 0) return(list(opinions = opinions, dt = Inf))  # no events possible

  # 2. Sample time to next event
  dt <- rexp(1, rate = R_tot)

  # 3. Select which node updates
  i <- sample(1:N, 1, prob = rates / R_tot)

  # 4. Node i adopts opinion of a random neighbor
  neighbors <- which(rig[i, ] != 0)

  if (length(neighbors) > 0) {

    j <- sample(neighbors, 1)
    opinions[i] <- opinions[j]

    }

  return(list(opinions = opinions, dt = dt))
}
