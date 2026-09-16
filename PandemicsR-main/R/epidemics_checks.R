#' Check Epidemic Position Indexes
#'
#' Validates the inverse position indexes used for O(1) removal from
#' the susceptible, infected, and recovered red/blue epidemic node lists.
#' For every node present in a list, its position index must point back
#' to the correct position in that list.
#'
#' @param state A simulation state object
#'
#' @return \code{TRUE} if all position indexes are valid.
#'
#' @examples
#' \dontrun{
#' check_epidemic_position_indexes(state)
#' }
check_epidemic_position_indexes <- function(state) {

  stopifnot(
    all(state$S_red_pos[state$S_red_nodes] ==
          seq_along(state$S_red_nodes)),

    all(state$S_blue_pos[state$S_blue_nodes] ==
          seq_along(state$S_blue_nodes)),

    all(state$I_red_pos[state$I_red_nodes] ==
          seq_along(state$I_red_nodes)),

    all(state$I_blue_pos[state$I_blue_nodes] ==
          seq_along(state$I_blue_nodes)),

    all(state$R_red_pos[state$R_red_nodes] ==
          seq_along(state$R_red_nodes)),

    all(state$R_blue_pos[state$R_blue_nodes] ==
          seq_along(state$R_blue_nodes)),

    all(state$S_pos[state$S_nodes] ==
          seq_along(state$S_nodes)),

    all(state$I_pos[state$I_nodes] ==
          seq_along(state$I_nodes)),

    all(state$R_pos[state$R_nodes] ==
          seq_along(state$R_nodes))
  )

  invisible(TRUE)
}
#' Check Epidemic Camp Membership
#'
#' Validates that recovered red and blue node lists are consistent with
#' the epidemic state and current opinions of individuals.
#'
#' For every recovered individual, this function checks that the individual
#' appears in the appropriate red or blue recovered node list according to
#' their current opinion.
#'
#' This function is intended primarily as a development and debugging
#' check for detecting inconsistencies between \code{state$epi},
#' \code{state$opinions}, \code{state$R_red_nodes}, and
#' \code{state$R_blue_nodes}.
#'
check_epidemic_camp_membership <- function(state) {

  if (!isTRUE(state$epidemic_started))
    return(invisible(TRUE))

  red_R <- state$R_red_nodes
  blue_R <- state$R_blue_nodes

  expected_red_R <- which(
    state$epi == state$R &
      state$opinions < 0
  )

  expected_blue_R <- which(
    state$epi == state$R &
      state$opinions > 0
  )

  if (!setequal(red_R, expected_red_R)) {

    cat("\n--- RED R MISMATCH ---\n")

    cat("R_red_nodes:",
        paste(red_R, collapse = ","), "\n")

    cat("expected_red_R:",
        paste(expected_red_R, collapse = ","), "\n")

    cat("In R_red_nodes but should not be:",
        paste(setdiff(red_R, expected_red_R), collapse = ","), "\n")

    cat("Should be R_red but missing:",
        paste(setdiff(expected_red_R, red_R), collapse = ","), "\n")
  }

  if (!setequal(blue_R, expected_blue_R)) {

    cat("\n--- BLUE R MISMATCH ---\n")

    cat("R_blue_nodes:",
        paste(blue_R, collapse = ","), "\n")

    cat("expected_blue_R:",
        paste(expected_blue_R, collapse = ","), "\n")

    cat("In R_blue_nodes but should not be:",
        paste(setdiff(blue_R, expected_blue_R), collapse = ","), "\n")

    cat("Should be R_blue but missing:",
        paste(setdiff(expected_blue_R, blue_R), collapse = ","), "\n")
  }

  invisible(TRUE)
}
