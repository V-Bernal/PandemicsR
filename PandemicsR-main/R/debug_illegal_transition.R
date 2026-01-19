#' Detect illegal Voter moves
#'
#' @docType package
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name bipartite_to_rig
#'
#' @import Matrix
#' @import igraph
#'
#' @references \url{}
#' @seealso \code{\link{brocolors}}
#' @keywords hplot
#'
#' @param old_opinions old_opinions
#' @param new_opinions new_opinions
#' @param t time
#'
#' @return NULL
#'
#' @examples
#' #---------------------------------
#' # debug_illegal_transition
#' #---------------------------------
#' @export
debug_illegal_transition <- function(old_opinions, new_opinions, t=NULL) {

  # allowed transitions:
  # extreme opinions (-2 and +2) MUST NOT change
  illegal <- which(
    (old_opinions %in% c(-2, 2)) &   # started extreme
      (new_opinions != old_opinions)   # changed to something else
  )

  if (length(illegal) > 0) {
    for (i in illegal) {
      cat("\n\n❌ ILLEGAL TRANSITION DETECTED at time =", t, "\n")
      cat("   Node:", i, "\n")
      cat("   Old:", old_opinions[i], "\n")
      cat("   New:", new_opinions[i], "\n")
      cat("   This should not happen (extreme opinions must stay fixed)\n")
      cat("------------------------------------------------------------\n")
    }
  }
}



