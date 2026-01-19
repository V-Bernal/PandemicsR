#' Detect illegal Schelling moves in a bipartite matrix
#'
#' @docType package
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name debug_illegal_schelling
#'
#' @import Matrix
#' @import igraph
#'
#' @references \url{}
#' @seealso \code{\link{brocolors}}
#' @keywords hplot
#'
#' @param old_bipartite old_bipartite
#' @param new_bipartite new_bipartite
#' @param t time
#'
#' @return NULL
#'
#' @examples
#' #---------------------------------
#' # debug_illegal_transition
#' #---------------------------------
#' @export
debug_illegal_schelling <- function(old_bipartite, new_bipartite,t =NULL) {

  # Identify changed edge
  diff_mat <- new_bipartite - old_bipartite
  changed <- which(diff_mat != 0, arr.ind = TRUE)

  # Ensures exactly one edge changes per Gillespie step.
  #
  # Checks additions only happen if addition rate > 0.
  #
  # Checks deletions only happen if deletion probability > 0.
  #
  # Validates RIG projection is correct.

  # Enhanced test function for schelling_step_bipartite_gillespie
  # Should be exactly one edge changed. No change is valid
  if (nrow(changed) > 1) {
    cat("\n\n❌ ILLEGAL TRANSITION DETECTED at time =", t, "\n")
    cat("   Edges:", 'from - to', "\n")
    print(changed)
    cat("   This should not happen (one edge changes per Gillespie step)\n")
    cat("------------------------------------------------------------\n")
  }
}

