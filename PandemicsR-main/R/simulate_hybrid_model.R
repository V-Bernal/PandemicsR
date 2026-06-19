#' Simulate Main-Compatible Voter-Schelling-SIR Model
#'
#' @author OpenAI Codex
#'
#' @name simulate_hybrid_model
#'
#' @param n Number of individuals.
#' @param m Number of groups.
#' @param t_max Maximum Gillespie time horizon.
#' @param lambda RIG weight parameter.
#' @param c_param Schelling group-joining rate parameter.
#' @param gamma Voter interaction rate used by the current main branch.
#' @param alpha Same-color radicalization rate.
#' @param alpha_deradicalization Same-color deradicalization rate.
#' @param alpha0 Spontaneous radicalization/deradicalization component.
#' @param beta_plus Schelling leave rate below threshold.
#' @param beta_minus Schelling leave rate above threshold.
#' @param T_threshold Minimum same-camp share to avoid the higher leave rate.
#' @param num_opinions Number of opinion states. Supports 2 and 4.
#' @param run_voter Logical toggle for voter dynamics.
#' @param run_schelling Logical toggle for Schelling dynamics.
#' @param run_epidemic Logical toggle for epidemic dynamics.
#' @param scaled_n Whether to divide the join rate by n.
#' @param scaled_m Whether to divide the join rate by m.
#' @param beta_red Infection rate for the red camp.
#' @param beta_blue Infection rate for the blue camp.
#' @param gamma_sir Recovery rate for the epidemic process.
#' @param initial_infected_fraction Initial infected fraction, or a count if
#'   greater than 1.
#' @param record_interval Gillespie time interval between recorded snapshots.
#' @param simulation_mode \code{"auto"}, \code{"exact"}, or \code{"aggregate"}.
#' @param exact_threshold Maximum \code{n} for exact individual simulation when
#'   \code{simulation_mode = "auto"}.
#' @param graph_threshold Maximum \code{n} for graph plotting metadata.
#' @param aggregate_max_steps Maximum integration steps for aggregate mode.
#'
#' @return A list containing graph objects, histories, and epidemic summaries.
#' @export
simulate_hybrid_model <- function(
    n, m, t_max, lambda, c_param,
    gamma = NULL,
    gamma_light = NULL,
    alpha = 0,
    alpha_deradicalization = 0,
    alpha0 = 0,
    beta_plus, beta_minus, T_threshold,
    num_opinions = 2,
    run_voter = TRUE, run_schelling = TRUE, run_epidemic = TRUE,
    scaled_n = FALSE, scaled_m = FALSE,
    beta_red = 0.5, beta_blue = 0.25, gamma_sir = 0.2,
    initial_infected_fraction = 0.05,
    record_every = NULL,
    record_interval = NULL,
    simulation_mode = c("auto", "exact", "aggregate"),
    exact_threshold = 500L,
    graph_threshold = 250L,
    aggregate_max_steps = 600L) {

  if (is.null(gamma)) {
    gamma <- gamma_light
  }

  params <- validate_hybrid_parameters(
    n = n,
    m = m,
    t_max = t_max,
    lambda = lambda,
    c_param = c_param,
    gamma = gamma,
    alpha = alpha,
    alpha_deradicalization = alpha_deradicalization,
    alpha0 = alpha0,
    beta_plus = beta_plus,
    beta_minus = beta_minus,
    T_threshold = T_threshold,
    num_opinions = num_opinions,
    beta_red = beta_red,
    beta_blue = beta_blue,
    gamma_sir = gamma_sir,
    initial_infected_fraction = initial_infected_fraction
  )
  list2env(params, environment())

  exact_threshold <- max(1L, as.integer(round(exact_threshold)))
  graph_threshold <- max(1L, as.integer(round(graph_threshold)))
  aggregate_max_steps <- max(1L, as.integer(round(aggregate_max_steps)))

  if (is.null(record_interval)) {
    record_interval <- if (!is.null(record_every)) as.numeric(record_every) else t_max / 500
  }
  record_interval <- max(.Machine$double.eps, as.numeric(record_interval))

  simulation_mode <- match.arg(simulation_mode)
  if (identical(simulation_mode, "aggregate") || (identical(simulation_mode, "auto") && n > exact_threshold)) {
    return(simulate_hybrid_aggregate_model(
      n = n,
      m = m,
      t_max = t_max,
      lambda = lambda,
      c_param = c_param,
      gamma = gamma,
      alpha = alpha,
      alpha_deradicalization = alpha_deradicalization,
      alpha0 = alpha0,
      beta_plus = beta_plus,
      beta_minus = beta_minus,
      T_threshold = T_threshold,
      num_opinions = num_opinions,
      run_voter = run_voter,
      run_schelling = run_schelling,
      run_epidemic = run_epidemic,
      scaled_n = scaled_n,
      scaled_m = scaled_m,
      beta_red = beta_red,
      beta_blue = beta_blue,
      gamma_sir = gamma_sir,
      initial_infected_fraction = initial_infected_fraction,
      max_steps = aggregate_max_steps
    ))
  }

  levels_vec <- get_levels_vec(num_opinions)
  state_labels <- get_state_labels(num_opinions)
  state_camps <- get_state_camp_index(num_opinions)
  n_states <- length(levels_vec)
  camp_labels <- get_camp_labels()
  camp_state_indices <- lapply(seq_along(camp_labels), function(idx) which(state_camps == idx))

  red_mod_idx <- match(-1, levels_vec)
  blue_mod_idx <- match(1, levels_vec)
  red_ext_idx <- match(-2, levels_vec)
  blue_ext_idx <- match(2, levels_vec)

  parse_initial_infected <- function(value, population_size) {
    if (!is.finite(value) || value <= 0) {
      return(0L)
    }
    if (value > 1) {
      return(min(population_size, as.integer(floor(value))))
    }
    min(population_size, as.integer(floor(value * population_size)))
  }

  ind_w <- rep(lambda, n)
  grp_w <- rep(lambda * n / m, m)

  B0 <- generate_bipartite(n, m, ind_w, grp_w)
  RIG0 <- bipartite_to_rig(B0)
  bipartite <- B0
  RIG <- RIG0
  rig_dirty <- FALSE

  opinions <- initialize_opinions_multi(n, num_opinions)
  opinion_history <- matrix(opinions, ncol = 1)

  groups_of_individual <- vector("list", n)
  for (i in seq_len(n)) {
    groups_of_individual[[i]] <- integer(0)
  }

  members <- vector("list", m)
  group_state_members <- vector("list", m)
  state_counts <- matrix(0L, nrow = m, ncol = n_states)

  for (g in seq_len(m)) {
    ids <- which(bipartite[, g] != 0)
    members[[g]] <- ids
    group_state_members[[g]] <- vector("list", n_states)
    for (state_idx in seq_len(n_states)) {
      group_state_members[[g]][[state_idx]] <- integer(0)
    }

    if (length(ids) > 0L) {
      state_idx <- match(opinions[ids], levels_vec)
      state_counts[g, ] <- tabulate(state_idx, nbins = n_states)

      for (state_id in seq_len(n_states)) {
        group_state_members[[g]][[state_id]] <- ids[state_idx == state_id]
      }

      for (i in ids) {
        groups_of_individual[[i]] <- c(groups_of_individual[[i]], g)
      }
    }
  }

  members0 <- lapply(members, identity)

  compartment_labels <- c("S", "I", "R")
  sir_state <- rep.int(1L, n)
  initial_infected_n <- if (isTRUE(run_epidemic)) parse_initial_infected(initial_infected_fraction, n) else 0L
  infected_seed <- integer(0L)
  if (initial_infected_n > 0L) {
    infected_seed <- if (initial_infected_n == n) seq_len(n) else sample.int(n, initial_infected_n)
    sir_state[infected_seed] <- 2L
  }

  compartment_state_members <- vector("list", length(compartment_labels))
  sir_counts <- matrix(
    0L,
    nrow = n_states,
    ncol = length(compartment_labels),
    dimnames = list(state_labels, compartment_labels)
  )

  for (comp_idx in seq_along(compartment_labels)) {
    compartment_state_members[[comp_idx]] <- vector("list", n_states)
    for (state_idx in seq_len(n_states)) {
      ids <- which(sir_state == comp_idx & opinions == levels_vec[[state_idx]])
      compartment_state_members[[comp_idx]][[state_idx]] <- ids
      sir_counts[state_idx, comp_idx] <- length(ids)
    }
  }

  sum_state_fraction <- function(current_opinions) {
    vapply(levels_vec, function(opinion_level) sum(current_opinions == opinion_level) / n, numeric(1))
  }

  sum_sir_counts <- function(current_counts) {
    c(
      S = sum(current_counts[, 1L]),
      I = sum(current_counts[, 2L]),
      R = sum(current_counts[, 3L])
    )
  }

  sum_camp_sir_counts <- function(current_counts) {
    camp_counts <- t(vapply(camp_state_indices, function(state_idx) {
      colSums(current_counts[state_idx, , drop = FALSE])
    }, numeric(length(compartment_labels))))
    rownames(camp_counts) <- camp_labels
    colnames(camp_counts) <- compartment_labels
    camp_counts
  }

  seed_camps <- camp_labels[state_camps[match(opinions[infected_seed], levels_vec)]]
  infection_time_history <- rep(0, length(infected_seed))
  infection_camp_history <- seed_camps

  record_infection <- function(chosen, current_time) {
    state_idx <- match(opinions[[chosen]], levels_vec)
    if (is.na(state_idx)) {
      return(invisible(FALSE))
    }
    infection_time_history <<- c(infection_time_history, current_time)
    infection_camp_history <<- c(infection_camp_history, camp_labels[[state_camps[[state_idx]]]])
    invisible(TRUE)
  }

  frac_history <- list(sum_state_fraction(opinions))
  sir_history <- list(sum_sir_counts(sir_counts))
  camp_sir_history <- list(sum_camp_sir_counts(sir_counts))
  time_history <- list(0)

  sample_weighted_index <- function(weights) {
    sample.int(length(weights), 1L, prob = weights / sum(weights))
  }

  sample_from_group_state <- function(group_idx, state_idx) {
    ids <- group_state_members[[group_idx]][[state_idx]]
    if (!length(ids)) {
      return(NA_integer_)
    }
    ids[[sample.int(length(ids), 1L)]]
  }

  sample_from_compartment_state <- function(comp_idx, state_idx) {
    ids <- compartment_state_members[[comp_idx]][[state_idx]]
    if (!length(ids)) {
      return(NA_integer_)
    }
    ids[[sample.int(length(ids), 1L)]]
  }

  move_opinion <- function(chosen, new_state_idx) {
    old_state_idx <- match(opinions[[chosen]], levels_vec)
    if (is.na(old_state_idx) || old_state_idx == new_state_idx) {
      return(invisible(FALSE))
    }

    current_compartment <- sir_state[[chosen]]
    opinions[[chosen]] <<- levels_vec[[new_state_idx]]

    compartment_state_members[[current_compartment]][[old_state_idx]] <<-
      remove_value_fast(compartment_state_members[[current_compartment]][[old_state_idx]], chosen)
    compartment_state_members[[current_compartment]][[new_state_idx]] <<-
      c(compartment_state_members[[current_compartment]][[new_state_idx]], chosen)

    sir_counts[old_state_idx, current_compartment] <<- sir_counts[old_state_idx, current_compartment] - 1L
    sir_counts[new_state_idx, current_compartment] <<- sir_counts[new_state_idx, current_compartment] + 1L

    for (g in groups_of_individual[[chosen]]) {
      group_state_members[[g]][[old_state_idx]] <<-
        remove_value_fast(group_state_members[[g]][[old_state_idx]], chosen)
      group_state_members[[g]][[new_state_idx]] <<-
        c(group_state_members[[g]][[new_state_idx]], chosen)
      state_counts[g, old_state_idx] <<- state_counts[g, old_state_idx] - 1L
      state_counts[g, new_state_idx] <<- state_counts[g, new_state_idx] + 1L
    }

    invisible(TRUE)
  }

  set_epidemic_state <- function(chosen, new_compartment) {
    old_compartment <- sir_state[[chosen]]
    if (old_compartment == new_compartment) {
      return(invisible(FALSE))
    }

    state_idx <- match(opinions[[chosen]], levels_vec)
    compartment_state_members[[old_compartment]][[state_idx]] <<-
      remove_value_fast(compartment_state_members[[old_compartment]][[state_idx]], chosen)
    compartment_state_members[[new_compartment]][[state_idx]] <<-
      c(compartment_state_members[[new_compartment]][[state_idx]], chosen)

    sir_counts[state_idx, old_compartment] <<- sir_counts[state_idx, old_compartment] - 1L
    sir_counts[state_idx, new_compartment] <<- sir_counts[state_idx, new_compartment] + 1L
    sir_state[[chosen]] <<- new_compartment
    invisible(TRUE)
  }

  add_membership <- function(chosen, group_idx) {
    state_idx <- match(opinions[[chosen]], levels_vec)
    members[[group_idx]] <<- c(members[[group_idx]], chosen)
    groups_of_individual[[chosen]] <<- c(groups_of_individual[[chosen]], group_idx)
    group_state_members[[group_idx]][[state_idx]] <<-
      c(group_state_members[[group_idx]][[state_idx]], chosen)
    state_counts[group_idx, state_idx] <<- state_counts[group_idx, state_idx] + 1L
    rig_dirty <<- TRUE
  }

  remove_membership <- function(chosen, group_idx) {
    state_idx <- match(opinions[[chosen]], levels_vec)
    members[[group_idx]] <<- remove_value_fast(members[[group_idx]], chosen)
    groups_of_individual[[chosen]] <<- remove_value_fast(groups_of_individual[[chosen]], group_idx)
    group_state_members[[group_idx]][[state_idx]] <<-
      remove_value_fast(group_state_members[[group_idx]][[state_idx]], chosen)
    state_counts[group_idx, state_idx] <<- state_counts[group_idx, state_idx] - 1L
    rig_dirty <<- TRUE
  }

  record_snapshot <- function(current_time) {
    frac_history[[length(frac_history) + 1L]] <<- sum_state_fraction(opinions)
    sir_history[[length(sir_history) + 1L]] <<- sum_sir_counts(sir_counts)
    camp_sir_history[[length(camp_sir_history) + 1L]] <<- sum_camp_sir_counts(sir_counts)
    time_history[[length(time_history) + 1L]] <<- current_time
  }

  t <- 0
  event_counter <- 0L
  last_record_time <- 0
  stop_reason <- "Reached maximum time"

  while (t < t_max) {
    if (isTRUE(run_epidemic) && sum(sir_counts[, 2L]) <= 0L) {
      stop_reason <- "No individual is infected (epidemic ended)"
      break
    }

    if (isTRUE(run_epidemic) && sum(sir_counts[, 3L]) == n) {
      stop_reason <- "All individuals recovered (epidemic ended)"
      break
    }

    if (num_opinions == 4L && all(opinions %in% c(-2, 2))) {
      stop_reason <- "Full polarization (all extreme opinions)"
      break
    }

    totals <- rowSums(state_counts)
    Ri <- if (!is.na(red_mod_idx)) state_counts[, red_mod_idx] else integer(m)
    Bi <- if (!is.na(blue_mod_idx)) state_counts[, blue_mod_idx] else integer(m)
    Ei_red <- if (!is.na(red_ext_idx)) state_counts[, red_ext_idx] else integer(m)
    Ei_blue <- if (!is.na(blue_ext_idx)) state_counts[, blue_ext_idx] else integer(m)
    Tot <- Ri + Bi

    enabled_moves <- integer(0)
    voter_term <- numeric(m)
    if (isTRUE(run_voter)) {
      enabled_moves <- c(enabled_moves, 1L, 2L)
      valid <- Ri > 0 & Bi > 0
      voter_term[valid] <- gamma * (Ri[valid] * Bi[valid] / Tot[valid])
    }

    join_term <- numeric(m)
    leaveR_rate <- numeric(m)
    leaveB_rate <- numeric(m)
    if (isTRUE(run_schelling)) {
      enabled_moves <- c(enabled_moves, 3L, 4L, 5L)
      join_term <- c_param * pmax(0, n - totals) * totals
      if (isTRUE(scaled_n)) {
        join_term <- join_term / n
      }
      if (isTRUE(scaled_m)) {
        join_term <- join_term / m
      }

      frac_red <- ifelse(Tot > 0, Ri / Tot, 0)
      frac_blue <- ifelse(Tot > 0, Bi / Tot, 0)
      leaveR_rate <- ifelse(frac_red < T_threshold, beta_plus * Ri, beta_minus * Ri)
      leaveB_rate <- ifelse(frac_blue < T_threshold, beta_plus * Bi, beta_minus * Bi)
    }

    infection_by_state <- numeric(n_states)
    recovery_by_state <- numeric(n_states)
    if (isTRUE(run_epidemic)) {
      infected_total <- sum(sir_counts[, 2L])
      susceptible_by_state <- sir_counts[, 1L]
      beta_by_state <- ifelse(state_camps == 1L, beta_red, beta_blue)
      infection_by_state <- if (infected_total > 0L) {
        infected_total / n * beta_by_state * susceptible_by_state
      } else {
        numeric(n_states)
      }
      recovery_by_state <- gamma_sir * sir_counts[, 2L]
    }
    infection_rate_tot <- sum(infection_by_state)
    recovery_rate_tot <- sum(recovery_by_state)

    Tot_g <- totals
    factor <- ifelse(Tot_g > 0, 1 / Tot_g, 0)
    rate_radicalize_red <- numeric(m)
    rate_radicalize_blue <- numeric(m)
    rate_deradicalize_red <- numeric(m)
    rate_deradicalize_blue <- numeric(m)
    if (num_opinions == 4L && isTRUE(run_voter)) {
      enabled_moves <- c(enabled_moves, 6L, 7L, 8L, 9L)
      rate_radicalize_red <- (alpha0 + alpha * Ei_red * factor) * Ri
      rate_radicalize_blue <- (alpha0 + alpha * Ei_blue * factor) * Bi
      rate_deradicalize_red <- (alpha0 + alpha_deradicalization * Ei_red * factor) * Ri
      rate_deradicalize_blue <- (alpha0 + alpha_deradicalization * Ei_blue * factor) * Bi
    }

    leaveR_rate_extreme <- numeric(m)
    leaveB_rate_extreme <- numeric(m)
    if (num_opinions == 4L && isTRUE(run_schelling)) {
      enabled_moves <- c(enabled_moves, 10L, 11L)
      frac_red_e <- ifelse(Tot_g > 0, Ei_red / Tot_g, 0)
      frac_blue_e <- ifelse(Tot_g > 0, Ei_blue / Tot_g, 0)
      leaveR_rate_extreme <- ifelse(frac_red_e < T_threshold, beta_plus * Ei_red, beta_minus * Ei_red)
      leaveB_rate_extreme <- ifelse(frac_blue_e < T_threshold, beta_plus * Ei_blue, beta_minus * Ei_blue)
    }

    lambda_i <- voter_term + join_term + leaveR_rate + leaveB_rate +
      rate_radicalize_red + rate_radicalize_blue +
      rate_deradicalize_red + rate_deradicalize_blue +
      leaveR_rate_extreme + leaveB_rate_extreme

    social_rate <- sum(lambda_i)
    epi_rate <- infection_rate_tot + recovery_rate_tot
    lambda_tot <- social_rate + epi_rate
    if (!is.finite(lambda_tot) || lambda_tot <= 0) {
      stop_reason <- "No enabled event has positive rate"
      break
    }

    t <- t + stats::rexp(1, lambda_tot)
    if (t >= t_max) {
      t <- t_max
      break
    }

    if (runif(1) < social_rate / lambda_tot) {
      if (social_rate <= 0) {
        next
      }

      group_idx <- sample_weighted_index(lambda_i)
      rates_vec <- c(
        ifelse(voter_term[group_idx] > 0, voter_term[group_idx] / 2, 0),
        ifelse(voter_term[group_idx] > 0, voter_term[group_idx] / 2, 0),
        join_term[group_idx],
        leaveR_rate[group_idx],
        leaveB_rate[group_idx],
        rate_radicalize_red[group_idx],
        rate_radicalize_blue[group_idx],
        rate_deradicalize_red[group_idx],
        rate_deradicalize_blue[group_idx],
        leaveR_rate_extreme[group_idx],
        leaveB_rate_extreme[group_idx]
      )
      rates_sub <- rates_vec[enabled_moves]
      if (!length(rates_sub) || sum(rates_sub) <= 0) {
        next
      }

      move <- sample(enabled_moves, 1L, prob = rates_sub / sum(rates_sub))
      if (move == 1L && Ri[[group_idx]] > 0 && Bi[[group_idx]] > 0) {
        chosen <- sample_from_group_state(group_idx, red_mod_idx)
        if (!is.na(chosen)) {
          move_opinion(chosen, blue_mod_idx)
        }
      } else if (move == 2L && Ri[[group_idx]] > 0 && Bi[[group_idx]] > 0) {
        chosen <- sample_from_group_state(group_idx, blue_mod_idx)
        if (!is.na(chosen)) {
          move_opinion(chosen, red_mod_idx)
        }
      } else if (move == 3L && totals[[group_idx]] < n) {
        chosen <- sample_outsider(group_idx, members, groups_of_individual, n)
        if (!is.na(chosen)) {
          add_membership(chosen, group_idx)
        }
      } else if (move == 4L && Ri[[group_idx]] > 0) {
        chosen <- sample_from_group_state(group_idx, red_mod_idx)
        if (!is.na(chosen)) {
          remove_membership(chosen, group_idx)
        }
      } else if (move == 5L && Bi[[group_idx]] > 0) {
        chosen <- sample_from_group_state(group_idx, blue_mod_idx)
        if (!is.na(chosen)) {
          remove_membership(chosen, group_idx)
        }
      } else if (move == 6L && Ei_red[[group_idx]] > 0 && Ri[[group_idx]] > 0) {
        chosen <- sample_from_group_state(group_idx, red_mod_idx)
        if (!is.na(chosen)) {
          move_opinion(chosen, red_ext_idx)
        }
      } else if (move == 7L && Ei_blue[[group_idx]] > 0 && Bi[[group_idx]] > 0) {
        chosen <- sample_from_group_state(group_idx, blue_mod_idx)
        if (!is.na(chosen)) {
          move_opinion(chosen, blue_ext_idx)
        }
      } else if (move == 8L && Ei_red[[group_idx]] > 0 && Ri[[group_idx]] > 0) {
        chosen <- sample_from_group_state(group_idx, red_ext_idx)
        if (!is.na(chosen)) {
          move_opinion(chosen, red_mod_idx)
        }
      } else if (move == 9L && Ei_blue[[group_idx]] > 0 && Bi[[group_idx]] > 0) {
        chosen <- sample_from_group_state(group_idx, blue_ext_idx)
        if (!is.na(chosen)) {
          move_opinion(chosen, blue_mod_idx)
        }
      } else if (move == 10L && Ei_red[[group_idx]] > 0) {
        chosen <- sample_from_group_state(group_idx, red_ext_idx)
        if (!is.na(chosen)) {
          remove_membership(chosen, group_idx)
        }
      } else if (move == 11L && Ei_blue[[group_idx]] > 0) {
        chosen <- sample_from_group_state(group_idx, blue_ext_idx)
        if (!is.na(chosen)) {
          remove_membership(chosen, group_idx)
        }
      }
    } else {
      if (epi_rate <= 0) {
        next
      }

      if (runif(1) < infection_rate_tot / epi_rate) {
        if (infection_rate_tot > 0) {
          state_idx <- sample_weighted_index(infection_by_state)
          chosen <- sample_from_compartment_state(1L, state_idx)
          if (!is.na(chosen) && isTRUE(set_epidemic_state(chosen, 2L))) {
            record_infection(chosen, t)
          }
        }
      } else if (recovery_rate_tot > 0) {
        state_idx <- sample_weighted_index(recovery_by_state)
        chosen <- sample_from_compartment_state(2L, state_idx)
        if (!is.na(chosen)) {
          set_epidemic_state(chosen, 3L)
        }
      }
    }

    event_counter <- event_counter + 1L
    if (t - last_record_time >= record_interval) {
      last_record_time <- t
      record_snapshot(t)
    }
  }

  opinion_history <- cbind(opinion_history, opinions)
  final_frac <- sum_state_fraction(opinions)
  final_sir <- sum_sir_counts(sir_counts)
  last_time <- time_history[[length(time_history)]]
  last_frac <- frac_history[[length(frac_history)]]
  last_sir <- sir_history[[length(sir_history)]]

  if (last_time < t || any(last_frac != final_frac) || any(last_sir != final_sir)) {
    record_snapshot(t)
  }

  frac_mat <- do.call(rbind, frac_history)
  colnames(frac_mat) <- state_labels

  sir_mat <- do.call(rbind, sir_history)
  colnames(sir_mat) <- compartment_labels

  time_history <- unlist(time_history, use.names = FALSE)
  infection_events <- data.frame(
    time = infection_time_history,
    camp = infection_camp_history,
    stringsAsFactors = FALSE
  )

  if (rig_dirty) {
    bipartite <- reconstruct_bipartite(members, n, m)
    RIG <- bipartite_to_rig(bipartite)
  }

  final_camp_sizes <- c(sum(opinions < 0), sum(opinions > 0))
  final_camp_susceptible <- c(
    sum(sir_state == 1L & opinions < 0),
    sum(sir_state == 1L & opinions > 0)
  )
  camp_attack_rate <- ifelse(final_camp_sizes > 0, 1 - final_camp_susceptible / final_camp_sizes, NA_real_)
  names(camp_attack_rate) <- camp_labels

  final_camp_sir <- rbind(
    c(
      S = sum(sir_state == 1L & opinions < 0),
      I = sum(sir_state == 2L & opinions < 0),
      R = sum(sir_state == 3L & opinions < 0)
    ),
    c(
      S = sum(sir_state == 1L & opinions > 0),
      I = sum(sir_state == 2L & opinions > 0),
      R = sum(sir_state == 3L & opinions > 0)
    )
  )
  rownames(final_camp_sir) <- camp_labels

  list(
    B0 = B0,
    B = bipartite,
    RIG0 = RIG0,
    RIG = RIG,
    opinions = opinions,
    opinion_history = opinion_history,
    num_opinions = num_opinions,
    members0 = members0,
    members = members,
    frac_mat = frac_mat,
    sir_mat = sir_mat,
    camp_sir_history = camp_sir_history,
    infection_events = infection_events,
    time_history = time_history,
    sir_state = sir_state,
    final_time = t,
    state_labels = state_labels,
    camp_labels = camp_labels,
    overall_attack_rate = 1 - sum(sir_state == 1L) / n,
    camp_attack_rate = camp_attack_rate,
    final_camp_sir = final_camp_sir,
    simulation_mode = "exact",
    graph_available = n <= graph_threshold,
    population_size = n,
    group_count = m,
    event_count = event_counter,
    stop_reason = stop_reason,
    model_notes = "Exact individual/network simulation using main-compatible event rates."
  )
}
