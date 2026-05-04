#' Build an alias table for O(1) discrete sampling
#'
#' Constructs the alias method data structure from a probability vector.
#' This allows sampling from a discrete distribution in constant time.
#'
#' @param p Numeric vector of probabilities (does not need to sum to 1;
#'   it will be normalized internally).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{prob}{Numeric vector of probabilities for each outcome}
#'   \item{alias}{Integer vector of alias indices}
#' }
#'
#' @details
#' The alias method preprocesses a discrete distribution in O(n) time
#' and allows subsequent sampling in O(1) time.
#'
#' @examples
#' #p <- c(0.1, 0.2, 0.7)
#' #alias <- build_alias(p)
#'
build_alias <- function(p) {
  n <- length(p)
  prob <- rep(0, n)
  alias <- rep(0, n)

  scaled <- p * n
  small <- which(scaled < 1)
  large <- which(scaled >= 1)

  while (length(small) && length(large)) {
    s <- small[1]; small <- small[-1]
    l <- large[1]; large <- large[-1]

    prob[s] <- scaled[s]
    alias[s] <- l

    scaled[l] <- scaled[l] - (1 - scaled[s])

    if (scaled[l] < 1) small <- c(small, l)
    else large <- c(large, l)
  }

  prob[large] <- 1
  prob[small] <- 1

  list(prob = prob, alias = alias)
}


#' Sample from an alias table in O(1) time
#'
#' Draws a single sample from a discrete distribution
#' represented by an alias table.
#'
#' @param alias_obj A list produced by \code{build_alias()},
#'   containing elements \code{prob} and \code{alias}.
#'
#' @return An integer index corresponding to the sampled outcome.
#'
#' @details
#' Each draw uses one uniform integer sample and one uniform
#' random number, achieving constant-time sampling.
#'
#' @examples
#' #p <- c(0.1, 0.2, 0.7)
#' #alias <- build_alias(p)
#' #alias_sample(alias)
#'
alias_sample <- function(alias_obj) {
  n <- length(alias_obj$prob)
  i <- sample.int(n, 1)
  if (runif(1) < alias_obj$prob[i]) i
  else alias_obj$alias[i]
}
