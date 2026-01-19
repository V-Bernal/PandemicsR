#' Single voter_step_gillespie_multi
#'
#' @docType package
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name voter_step_gillespie_multi
#'
#' @import Matrix
#' @import igraph
#'
#' @param rig number of individuals vertices
#' @param opinions number of group vertices
#' @param kappa Poisson clock ringing rate (see \code{\link[graphics]{par}})
#' @param decision_matrix decision_matrix
#' @param vertex vertex
#' @return opinions, dt \code{opinions}.
#'
#' @examples
#' #---------------------------------
#' # voter_step_gillespie_multi
#' # voter_step_gillespie_multi (rig, opinions, kappa)
#' #---------------------------------
#' @export
voter_step_gillespie_multi <- function(rig, opinions, vertex = NULL, kappa, decision_matrix =NULL){
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

  # if (length(neighbors) > 0) {
  #
  #   j <- sample(neighbors, 1)
  #   opinions[i] <- opinions[j]
  #
  #   }
  if (length(neighbors) > 1){

    j <- sample(x = neighbors, size = 1, replace = FALSE)

  } else if(length(neighbors) == 1) { # Watch out*

    j <- neighbors
  } else{ return(opinions)}

  #========================
  # 4. Update opinion according to labels decision matrix
  levels_vec <- get_levels_vec(ncol(decision_matrix))#c(-2, -1, 1, 2)
  id <- as.integer(factor(opinions, levels = levels_vec))
  current_op <- id[i]
  neighbor_op <- id[j]
  prob <- decision_matrix[current_op, neighbor_op]
  #========================

  if (is.na(current_op) | is.na(neighbor_op)) {
    stop("Illegal opinion value detected: ", opinions[i], " or ", opinions[j])
  }


  # 4. Update opinion according to index i, j in the decision matrix
  #prob <- decision_matrix[ i, j]
  if (runif(1) < prob) {

    #debug_illegal_transition(opinions[i], opinions[j], t = NULL)

    opinions[i] <- opinions[j]# neighbor_op
  }

  return(list(opinions = opinions, dt = dt))
}
