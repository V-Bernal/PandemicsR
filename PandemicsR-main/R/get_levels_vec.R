#' Generate Opinion Levels
#'
#' Creates a vector of discrete opinion levels based on the specified number
#' of opinions. The levels are evenly spaced around 0
#'
#' @author Victor Bernal
#'
#' @name get_levels_vec
#'
#' @param num_opinions Number of opinions
#'
#' @return A numeric vector of length \code{num_opinions} around 0
#'
#' @details
#' The levels are generated using an evenly spaced sequence
#'
#' @examples
#' # Generate 3 opinion levels
#' get_levels_vec(3)
#' # Returns: -2 -1  1  2
#'
#' @export
get_levels_vec <- function(num_opinions){

  lv <- seq(-num_opinions/2, num_opinions/2, by = 1)

  lv <- setdiff(lv, 0)

  return(lv)
}
