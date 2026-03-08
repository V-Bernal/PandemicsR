#' make_transition_matrix
#'
#' @docType package
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name make_transition_matrix
#'
#' @import Matrix
#' @import igraph
#'
#' @references \url{}
#' @seealso \code{\link{brocolors}}
#' @keywords hplot
#'
#' @param N
#'
#' @return transition matrix
#'
#' @examples
#' #---------------------------------
#' # make_transition_matrix
#' #---------------------------------
#' @export
make_transition_matrix <- function(N) {

  M <- matrix(0, nrow = N, ncol = N)

  for (i in 2:(N-1)) {
    # distances from opinion i
    d <- abs(i - 1:N)

    # weights (higher for closer opinions)
    w <- (N - d)

    # no self-transition
    w[i] <- 0

    # normalize to sum to 1
    M[i, ] <- w / sum(w)
  }

  # extremes M[1, ] and M[N, ] remain all zeros
  if( !any(rowSums(M) %in% c(0,1)) ){
    cat('Invalid transition matrix')
    break
  }

     return(M)
}
