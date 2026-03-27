#' Initialize Voters' Opinions (Multiple States)
#'
#' Generates an initial configuration of opinions for a population of
#' individuals, where each individual is assigned one of several possible
#' discrete opinion states.
#'
#' @author Victor Bernal
#'
#' @name initialize_opinions_multi
#'
#' @import Matrix
#' @import igraph
#'
#' @param n Integer. Number of individual vertices.
#' @param num_opinions Integer. Number of distinct opinion states.
#'
#' @return An integer vector of length \code{n}, where each entry represents
#' the opinion of an individual.
#'
#' @details
#' Each individual is assigned an opinion by sampling independently from a
#' discrete uniform distribution over \code{1:num_opinions}. This produces
#' an unbiased initial distribution of opinions across the population.
#'
#' @examples
#' #---------------------------------
#' # Initialize opinions
#' # initialize_opinions_multi(n = 15, num_Opinions = 2)
#' #---------------------------------
#' @export
initialize_opinions_multi <- function(n, num_opinions = NULL) {

  set_opinions <- get_levels_vec(num_opinions)

  opinions_sample <- sample(set_opinions , n, replace = TRUE)

  return(opinions_sample)

  }
