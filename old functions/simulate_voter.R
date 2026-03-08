#' Simulation voter's dynamics
#'
#' @docType package
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name simulate_voter
#'
#' @import Matrix
#' @import igraph
#'
#' @param rig number of individuals vertices
#' @param opinions number of group vertices
#' @param times number of iterations
#' @param frac_pos fraction of positive opinions +1
#' @param kappa Poisson clock ringing rate (see \code{\link[graphics]{par}})
#'
#' @return opinions \code{opinions}.
#'
#' @examples
#' #---------------------------------
#' # Simulate Voter
#' # simulate_voter(rig, opinions, times, kappa)
#' #---------------------------------
#' @export
#'
simulate_voter <- function(rig, opinions, times, kappa) {
  # This function
  frac_pos <- as.vector(NA*numeric(times))  # fraction of +1 opinions per iteration
  opinion_history <- matrix(NA, nrow = length(opinions), ncol = times)
  rownames(opinion_history) <- 1:length(opinions)

  for (t in 1:times){

    # opinions_temp <- voter_step(rig, opinions, kappa, vertex)
    # #opinions <- opinions_temp[[1]]
    # #opinion_history <- opinions_temp[[2]]
    # #frac_pos <- opinions_temp[[3]]

    # Update opinions
    opinions <- voter_step(rig, opinions, kappa, vertex)

    # Record history
    opinion_history[, t] <- opinions

    # Compute fraction of +1 opinions
    frac_pos[t] <- mean(opinions == 1)

  }

  return(list(opinions = opinions,
                     opinion_history = opinion_history,
                     frac_pos = frac_pos))
}
