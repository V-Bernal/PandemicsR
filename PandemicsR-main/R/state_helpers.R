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

#' Validate Hybrid Model Parameters
#'
#' @author OpenAI Codex
#'
#' @name validate_hybrid_parameters
#'
#' @param n Number of individuals.
#' @param m Number of groups.
#' @param t_max Maximum Gillespie time horizon.
#' @param lambda RIG weight parameter.
#' @param c_param Schelling group-joining rate parameter.
#' @param gamma_light Baseline voter rate for light opinion states.
#' @param gamma_dark Baseline voter rate for dark opinion states.
#' @param infected_dark_multiplier Multiplier applied to dark-state voter rates
#'   when the individual is infected.
#' @param beta_plus Schelling leave rate below threshold.
#' @param beta_minus Schelling leave rate above threshold.
#' @param T_threshold Minimum same-camp share to avoid the higher leave rate.
#' @param num_opinions Number of opinion states. Supports 2 and 4.
#' @param beta_red Infection rate for the red camp.
#' @param beta_blue Infection rate for the blue camp.
#' @param gamma_sir Recovery rate for the epidemic process.
#' @param initial_infected_fraction Initial infected fraction, or a count if
#'   greater than 1.
#'
#' @return A normalized list of scalar model parameters.
#' @export
validate_hybrid_parameters <- function(
    n, m, t_max, lambda, c_param,
    gamma_light, gamma_dark,
    infected_dark_multiplier,
    beta_plus, beta_minus, T_threshold,
    num_opinions,
    beta_red, beta_blue, gamma_sir,
    initial_infected_fraction) {

  require_scalar <- function(value, name) {
    if (length(value) != 1L || !is.finite(value)) {
      stop(sprintf("%s must be a finite scalar.", name), call. = FALSE)
    }
    as.numeric(value)
  }

  require_nonnegative <- function(value, name) {
    value <- require_scalar(value, name)
    if (value < 0) {
      stop(sprintf("%s must be non-negative.", name), call. = FALSE)
    }
    value
  }

  require_probability <- function(value, name) {
    value <- require_scalar(value, name)
    if (value < 0 || value > 1) {
      stop(sprintf("%s must be between 0 and 1.", name), call. = FALSE)
    }
    value
  }

  require_count <- function(value, name, min_value) {
    value <- require_scalar(value, name)
    rounded <- as.integer(round(value))
    if (abs(value - rounded) > .Machine$double.eps^0.5 || rounded < min_value) {
      stop(sprintf("%s must be an integer greater than or equal to %s.", name, min_value), call. = FALSE)
    }
    rounded
  }

  n <- require_count(n, "n", 1L)
  m <- require_count(m, "m", 1L)
  if (m > n) {
    stop("m cannot be larger than n.", call. = FALSE)
  }

  num_opinions <- require_count(num_opinions, "num_opinions", 2L)
  if (!(num_opinions %in% c(2L, 4L))) {
    stop("This implementation currently supports 2 or 4 opinion states.", call. = FALSE)
  }

  list(
    n = n,
    m = m,
    t_max = require_nonnegative(t_max, "t_max"),
    lambda = require_nonnegative(lambda, "lambda"),
    c_param = require_probability(c_param, "c_param"),
    gamma_light = require_nonnegative(gamma_light, "gamma_light"),
    gamma_dark = require_nonnegative(gamma_dark, "gamma_dark"),
    infected_dark_multiplier = require_nonnegative(infected_dark_multiplier, "infected_dark_multiplier"),
    beta_plus = require_probability(beta_plus, "beta_plus"),
    beta_minus = require_probability(beta_minus, "beta_minus"),
    T_threshold = require_probability(T_threshold, "T_threshold"),
    num_opinions = num_opinions,
    beta_red = require_nonnegative(beta_red, "beta_red"),
    beta_blue = require_nonnegative(beta_blue, "beta_blue"),
    gamma_sir = require_nonnegative(gamma_sir, "gamma_sir"),
    initial_infected_fraction = require_nonnegative(initial_infected_fraction, "initial_infected_fraction")
  )
}
