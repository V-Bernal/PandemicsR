#' Single Schelling step on the bipartite graph with gillespie
#'
#' @docType package
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name schelling_step_bipartite_gillespie_multi
#'
#' @import Matrix
#' @import igraph
#'
#' @param bipartite_matrix bipartite grpah
#' @param opinions number of group vertices
#' @param RIG RIG
#' @param const cont
#'
#' @return bipartite_matrix, RIG, dt
#'
#' @examples
#' #---------------------------------
#' # A single Schelling step
#' #schelling_step_bipartite_gillespie_multi (bipartite_matrix, RIG, opinions, const = 1)
#' #---------------------------------
#' @export
schelling_step_bipartite_gillespie_multi <- function(bipartite_matrix, RIG, opinions, const = 1) {

  # Gillespie updates exactly ONE event per step:
  # Event = add or remove one bipartite edge (i,j)
  # with probability proportional to its event rate.
  #
  # So the new Gillespie step must:
  # 1 List all possible events
  # 2 Assign an event rate to each
  # 3 Compute total rate R
  # 4 Sample the time to next event
  # 5 Choose exactly ONE event proportional to its rate
  # 6 Execute that event (add/remove)
  # 7 Update the bipartite matrix and RIG
  # 8 Return the new configuration + dt


  n <- nrow(bipartite_matrix)
  m <- ncol(bipartite_matrix)
  Nnodes <- n

  p_addition <- const / m  # constant rate

  # 1 List all possible events. Assign an event rate to each
  event_list <- list()
  rates <- c()

  idx <- 1

  # Each pair (i,j) can create exactly one event:
    # Loop over all node–object pairs (i,j)
    for (i in 1:n) {
      alpha <- opinions[i]

      for (j in 1:m) {

        # Count positive / negative neighbors at object j
        # k_pos <- sum(bipartite_matrix[, j] * (opinions == 1))
        # k_neg <- sum(bipartite_matrix[, j] * (opinions == -1))
        k_pos <- sum(bipartite_matrix[, j] * (opinions > 0))
        k_neg <- sum(bipartite_matrix[, j] * (opinions < 0 ))


        # Avoid division by zero
        if ((k_pos + k_neg) == 0) next

        if (bipartite_matrix[i, j] == 0) {
          # ADDITION event
          rate <- p_addition

          event_list[[idx]] <- list(type="add", i=i, j=j)
          rates[idx] <- rate
          idx <- idx + 1

        } else {
          # DELETION event
          p_del <- -1 * alpha * (k_pos - k_neg) / (k_pos + k_neg)
          p_del <- max(p_del, 0)   # enforce non-negative rate

          event_list[[idx]] <- list(type="del", i=i, j=j)
          rates[idx] <- p_del
          idx <- idx + 1
        }
      }
    }

  # Compute total rate R
  R_tot <- sum(rates)

  #error handling
  if (R_tot == 0) {
    return(list(
      bipartite_matrix = bipartite_matrix,
      RIG = RIG,
      dt = Inf
    ))
  }

  # 2. Sample waiting time (time to next event)
  dt <- rexp(1, rate = R_tot)

  # 3. Choose exactly ONE event proportional to its rate
  chosen <- sample(1:length(rates), 1, prob = rates / R_tot)
  ev <- event_list[[chosen]]

  # 4. Execute that event (add/remove)
  if (ev$type == "add") {
    bipartite_matrix[ev$i, ev$j] <- 1
  } else {
    bipartite_matrix[ev$i, ev$j] <- 0
  }

  # 5. Update the bipartite matrix and RIG
  RIG <- bipartite_to_rig(bipartite_matrix)

  return(list(
    bipartite_matrix = bipartite_matrix,
    RIG = RIG,
    dt = dt
  ))
}
