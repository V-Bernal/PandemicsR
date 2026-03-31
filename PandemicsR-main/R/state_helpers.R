#' State Labels
#'
#' @author OpenAI Codex
#'
#' @name get_state_labels
#'
#' @param num_opinions Number of opinion states.
#'
#' @return Character vector of labels aligned with \code{get_levels_vec()}.
#' @export
get_state_labels <- function(num_opinions) {
  levels_vec <- get_levels_vec(num_opinions)

  if (identical(levels_vec, c(-2, -1, 1, 2))) {
    return(c("Dark red", "Light red", "Light blue", "Dark blue"))
  }

  if (identical(levels_vec, c(-1, 1))) {
    return(c("Red", "Blue"))
  }

  as.character(levels_vec)
}

#' State Camp Index
#'
#' @author OpenAI Codex
#'
#' @name get_state_camp_index
#'
#' @param num_opinions Number of opinion states.
#'
#' @return Integer vector where 1 = red camp and 2 = blue camp.
#' @export
get_state_camp_index <- function(num_opinions) {
  ifelse(get_levels_vec(num_opinions) < 0, 1L, 2L)
}

#' Camp Labels
#'
#' @author OpenAI Codex
#'
#' @name get_camp_labels
#'
#' @return Character vector with red and blue camp labels.
#' @export
get_camp_labels <- function() {
  c("Red camp", "Blue camp")
}

#' Dark-State Flags
#'
#' @author OpenAI Codex
#'
#' @name get_dark_state_flags
#'
#' @param num_opinions Number of opinion states.
#'
#' @return Logical vector indicating which states are dark / stubborn.
#' @export
get_dark_state_flags <- function(num_opinions) {
  levels_vec <- get_levels_vec(num_opinions)
  if (length(levels_vec) <= 2L) {
    return(rep(FALSE, length(levels_vec)))
  }

  abs(levels_vec) == max(abs(levels_vec))
}

#' Opinion Step
#'
#' @author OpenAI Codex
#'
#' @name opinion_step
#'
#' @param opinion Current opinion value.
#' @param direction Integer step direction, \code{-1} or \code{+1}.
#' @param num_opinions Number of opinion states.
#'
#' @return Updated opinion value after one adjacent move.
#' @export
opinion_step <- function(opinion, direction, num_opinions) {
  levels_vec <- get_levels_vec(num_opinions)
  idx <- match(opinion, levels_vec)

  if (is.na(idx)) {
    stop("Unknown opinion value: ", opinion)
  }

  new_idx <- min(length(levels_vec), max(1L, idx + direction))
  levels_vec[[new_idx]]
}

#' Opinion State Index
#'
#' @author OpenAI Codex
#'
#' @name opinion_state_index
#'
#' @param opinions Opinion values.
#' @param num_opinions Number of opinion states.
#'
#' @return Integer vector of state indices.
#' @export
opinion_state_index <- function(opinions, num_opinions) {
  match(opinions, get_levels_vec(num_opinions))
}

#' Remove a Value From a Vector
#'
#' @author OpenAI Codex
#'
#' @name remove_value_fast
#'
#' @param vec Integer vector.
#' @param value Value to remove.
#'
#' @return Input vector with the first match removed.
#' @export
remove_value_fast <- function(vec, value) {
  idx <- match(value, vec, nomatch = 0L)
  if (idx == 0L) {
    return(vec)
  }

  last_idx <- length(vec)
  vec[idx] <- vec[last_idx]
  vec[-last_idx]
}

#' Reconstruct a Bipartite Matrix From Membership Lists
#'
#' @author OpenAI Codex
#'
#' @name reconstruct_bipartite
#'
#' @param members List of group memberships.
#' @param n Number of individuals.
#' @param m Number of groups.
#'
#' @return Sparse bipartite matrix.
#' @export
reconstruct_bipartite <- function(members, n, m) {
  group_sizes <- lengths(members)
  if (!any(group_sizes)) {
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0), x = integer(0), dims = c(n, m)))
  }

  Matrix::sparseMatrix(
    i = unlist(members, use.names = FALSE),
    j = rep.int(seq_len(m), group_sizes),
    x = 1L,
    dims = c(n, m)
  )
}

#' Sample an Outsider For a Group
#'
#' @author OpenAI Codex
#'
#' @name sample_outsider
#'
#' @param group_i Group index.
#' @param members Membership list by group.
#' @param groups_of_individual Reverse membership index.
#' @param n Number of individuals.
#' @param max_tries Number of rejection-sampling attempts.
#'
#' @return Integer individual id or \code{NA}.
#' @export
sample_outsider <- function(group_i, members, groups_of_individual, n, max_tries = 25L) {
  for (attempt in seq_len(max_tries)) {
    candidate <- sample.int(n, 1L)
    if (!(group_i %in% groups_of_individual[[candidate]])) {
      return(candidate)
    }
  }

  available <- setdiff(seq_len(n), members[[group_i]])
  if (!length(available)) {
    return(NA_integer_)
  }

  if (length(available) == 1L) {
    return(available)
  }

  sample(available, 1L)
}
