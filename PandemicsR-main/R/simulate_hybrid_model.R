#' Simulate Hybrid Voter-Schelling-SIR Model
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
#' @param gamma_light Baseline voter rate for light opinion states.
#' @param gamma_dark Baseline voter rate for dark opinion states.
#' @param infected_dark_multiplier Multiplier applied to dark-state voter
#'   rates when the individual is infected.
#' @param beta_plus Schelling leave rate below threshold.
#' @param beta_minus Schelling leave rate above threshold.
#' @param T_threshold Minimum same-camp share to avoid the higher leave rate.
#' @param num_opinions Number of opinion states. Supports 2 and 4.
#' @param run_voter Logical toggle for voter dynamics.
#' @param run_schelling Logical toggle for Schelling dynamics.
#' @param beta_red Infection rate for the red camp.
#' @param beta_blue Infection rate for the blue camp.
#' @param gamma_sir Recovery rate for the epidemic process.
#' @param initial_infected_fraction Initial infected fraction, or a count if
#'   greater than 1.
#' @param record_every Number of Gillespie events between recorded snapshots.
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
    gamma_light, gamma_dark = gamma_light,
    infected_dark_multiplier = 1,
    beta_plus, beta_minus, T_threshold,
    num_opinions = 2,
    run_voter = TRUE, run_schelling = TRUE,
    beta_red = 0.5, beta_blue = 0.25, gamma_sir = 0.2,
    initial_infected_fraction = 0.05,
    record_every = max(1L, as.integer(n)),
    simulation_mode = c("auto", "exact", "aggregate"),
    exact_threshold = 500L,
    graph_threshold = 250L,
    aggregate_max_steps = 600L) {

  params <- validate_hybrid_parameters(
    n = n,
    m = m,
    t_max = t_max,
    lambda = lambda,
    c_param = c_param,
    gamma_light = gamma_light,
    gamma_dark = gamma_dark,
    infected_dark_multiplier = infected_dark_multiplier,
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
  record_every <- max(1L, as.integer(round(record_every)))
  exact_threshold <- max(1L, as.integer(round(exact_threshold)))
  graph_threshold <- max(1L, as.integer(round(graph_threshold)))
  aggregate_max_steps <- max(1L, as.integer(round(aggregate_max_steps)))

  simulation_mode <- match.arg(simulation_mode)
  if (identical(simulation_mode, "aggregate") || (identical(simulation_mode, "auto") && n > exact_threshold)) {
    return(simulate_hybrid_aggregate_model(
      n = n,
      m = m,
      t_max = t_max,
      lambda = lambda,
      c_param = c_param,
      gamma_light = gamma_light,
      gamma_dark = gamma_dark,
      infected_dark_multiplier = infected_dark_multiplier,
      beta_plus = beta_plus,
      beta_minus = beta_minus,
      T_threshold = T_threshold,
      num_opinions = num_opinions,
      run_voter = run_voter,
      run_schelling = run_schelling,
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
  dark_states <- get_dark_state_flags(num_opinions)
  n_states <- length(levels_vec)
  camp_labels <- get_camp_labels()
  camp_state_indices <- lapply(seq_along(camp_labels), function(idx) which(state_camps == idx))

  if (num_opinions <= 2L) {
    gamma_dark <- gamma_light
    infected_dark_multiplier <- 1
  }

  base_state_rates <- rep(gamma_light, n_states)
  base_state_rates[dark_states] <- gamma_dark
  infected_state_rates <- base_state_rates
  infected_state_rates[dark_states] <- infected_state_rates[dark_states] * infected_dark_multiplier

  ind_w <- rep(lambda, n)
  grp_w <- rep(lambda * n / m, m)

  B0 <- generate_bipartite(n, m, ind_w, grp_w)
  RIG0 <- bipartite_to_rig(B0)
  bipartite <- B0
  RIG <- RIG0
  rig_dirty <- FALSE

  opinions <- initialize_opinions_multi(n, num_opinions)
  opinion_history <- matrix(opinions, ncol = 1)

  parse_initial_infected <- function(value, population_size) {
    if (!is.finite(value) || value <= 0) {
      return(0L)
    }

    if (value > 1) {
      return(min(population_size, as.integer(round(value))))
    }

    min(population_size, max(1L, as.integer(round(value * population_size))))
  }

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

  sir_state <- rep.int(1L, n)
  initial_infected_n <- parse_initial_infected(initial_infected_fraction, n)
  infected_seed <- integer(0L)
  if (initial_infected_n > 0L) {
    infected_seed <- if (initial_infected_n == n) seq_len(n) else sample.int(n, initial_infected_n)
    sir_state[infected_seed] <- 2L
  }

  compartment_labels <- c("S", "I", "R")
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

  infected_counts <- matrix(0L, nrow = m, ncol = n_states)
  for (g in seq_len(m)) {
    ids <- members[[g]]
    if (!length(ids)) {
      next
    }

    infected_ids <- ids[sir_state[ids] == 2L]
    if (length(infected_ids)) {
      infected_state_idx <- match(opinions[infected_ids], levels_vec)
      infected_counts[g, ] <- tabulate(infected_state_idx, nbins = n_states)
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

  sample_weighted_state <- function(weights) {
    idx <- sample.int(length(weights), 1L, prob = weights / sum(weights))
    idx
  }

  sample_from_group_state <- function(group_idx, state_idx) {
    ids <- group_state_members[[group_idx]][[state_idx]]
    if (!length(ids)) {
      return(NA_integer_)
    }

    if (length(ids) == 1L || !dark_states[[state_idx]]) {
      return(ids[[sample.int(length(ids), 1L)]])
    }

    weights <- ifelse(sir_state[ids] == 2L, infected_dark_multiplier, 1)
    if (!all(is.finite(weights)) || sum(weights) <= 0) {
      return(ids[[sample.int(length(ids), 1L)]])
    }

    sample(ids, 1L, prob = weights)
  }

  sample_from_compartment_state <- function(comp_idx, state_idx) {
    ids <- compartment_state_members[[comp_idx]][[state_idx]]
    if (!length(ids)) {
      return(NA_integer_)
    }

    if (length(ids) == 1L) {
      return(ids)
    }

    sample(ids, 1L)
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

      if (current_compartment == 2L) {
        infected_counts[g, old_state_idx] <<- infected_counts[g, old_state_idx] - 1L
        infected_counts[g, new_state_idx] <<- infected_counts[g, new_state_idx] + 1L
      }
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

    if (old_compartment == 2L || new_compartment == 2L) {
      delta <- if (new_compartment == 2L) 1L else -1L
      for (g in groups_of_individual[[chosen]]) {
        infected_counts[g, state_idx] <<- infected_counts[g, state_idx] + delta
      }
    }

    invisible(TRUE)
  }

  add_membership <- function(chosen, group_idx) {
    state_idx <- match(opinions[[chosen]], levels_vec)
    members[[group_idx]] <<- c(members[[group_idx]], chosen)
    groups_of_individual[[chosen]] <<- c(groups_of_individual[[chosen]], group_idx)
    group_state_members[[group_idx]][[state_idx]] <<-
      c(group_state_members[[group_idx]][[state_idx]], chosen)
    state_counts[group_idx, state_idx] <<- state_counts[group_idx, state_idx] + 1L
    if (sir_state[[chosen]] == 2L) {
      infected_counts[group_idx, state_idx] <<- infected_counts[group_idx, state_idx] + 1L
    }
    rig_dirty <<- TRUE
  }

  remove_membership <- function(chosen, group_idx) {
    state_idx <- match(opinions[[chosen]], levels_vec)
    members[[group_idx]] <<- remove_value_fast(members[[group_idx]], chosen)
    groups_of_individual[[chosen]] <<- remove_value_fast(groups_of_individual[[chosen]], group_idx)
    group_state_members[[group_idx]][[state_idx]] <<-
      remove_value_fast(group_state_members[[group_idx]][[state_idx]], chosen)
    state_counts[group_idx, state_idx] <<- state_counts[group_idx, state_idx] - 1L
    if (sir_state[[chosen]] == 2L) {
      infected_counts[group_idx, state_idx] <<- infected_counts[group_idx, state_idx] - 1L
    }
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

  while (t < t_max) {
    totals <- rowSums(state_counts)
    red_counts <- rowSums(state_counts[, camp_state_indices[[1]], drop = FALSE])
    blue_counts <- rowSums(state_counts[, camp_state_indices[[2]], drop = FALSE])

    voter_down <- matrix(0, nrow = m, ncol = n_states)
    voter_up <- matrix(0, nrow = m, ncol = n_states)

    if (run_voter) {
      for (g in seq_len(m)) {
        total_g <- totals[[g]]
        if (total_g < 2L) {
          next
        }

        counts_g <- state_counts[g, ]
        cumulative <- cumsum(counts_g)
        lower_counts <- c(0L, cumulative[-length(cumulative)])
        higher_counts <- total_g - cumulative

        for (state_idx in which(counts_g > 0L)) {
          noninfected_count <- counts_g[[state_idx]] - infected_counts[g, state_idx]
          weighted_changers <- base_state_rates[[state_idx]] * noninfected_count +
            infected_state_rates[[state_idx]] * infected_counts[g, state_idx]

          if (lower_counts[[state_idx]] > 0L) {
            voter_down[g, state_idx] <- weighted_changers * lower_counts[[state_idx]] / total_g
          }

          if (higher_counts[[state_idx]] > 0L) {
            voter_up[g, state_idx] <- weighted_changers * higher_counts[[state_idx]] / total_g
          }
        }
      }
    }

    join_term <- numeric(m)
    leave_red_rate <- numeric(m)
    leave_blue_rate <- numeric(m)

    if (run_schelling) {
      join_term <- (c_param / m) * pmax(0, n - totals)
      red_frac <- ifelse(totals > 0L, red_counts / totals, 0)
      blue_frac <- ifelse(totals > 0L, blue_counts / totals, 0)
      leave_red_rate <- ifelse(red_frac < T_threshold, beta_plus * red_counts, beta_minus * red_counts)
      leave_blue_rate <- ifelse(blue_frac < T_threshold, beta_plus * blue_counts, beta_minus * blue_counts)
    }

    infected_total <- sum(sir_counts[, 2L])
    susceptible_red <- sum(sir_counts[camp_state_indices[[1]], 1L])
    susceptible_blue <- sum(sir_counts[camp_state_indices[[2]], 1L])
    infection_red <- if (infected_total > 0L) beta_red * infected_total / n * susceptible_red else 0
    infection_blue <- if (infected_total > 0L) beta_blue * infected_total / n * susceptible_blue else 0
    recovery_by_state <- gamma_sir * sir_counts[, 2L]

    class_rates <- c(
      voter_down = sum(voter_down),
      voter_up = sum(voter_up),
      join = sum(join_term),
      leave_red = sum(leave_red_rate),
      leave_blue = sum(leave_blue_rate),
      infect_red = infection_red,
      infect_blue = infection_blue,
      recover = sum(recovery_by_state)
    )

    lambda_tot <- sum(class_rates)
    if (!is.finite(lambda_tot) || lambda_tot <= 0) {
      break
    }

    dt <- stats::rexp(1, lambda_tot)
    t <- t + dt
    if (t >= t_max) {
      t <- t_max
      break
    }

    event_type <- names(class_rates)[sample.int(length(class_rates), 1L, prob = class_rates / lambda_tot)]

    if (identical(event_type, "voter_down")) {
      weights <- as.vector(voter_down)
      picked <- arrayInd(sample.int(length(weights), 1L, prob = weights / sum(weights)), .dim = dim(voter_down))
      group_idx <- picked[[1]]
      state_idx <- picked[[2]]
      chosen <- sample_from_group_state(group_idx, state_idx)
      if (!is.na(chosen)) {
        move_opinion(chosen, state_idx - 1L)
      }
    } else if (identical(event_type, "voter_up")) {
      weights <- as.vector(voter_up)
      picked <- arrayInd(sample.int(length(weights), 1L, prob = weights / sum(weights)), .dim = dim(voter_up))
      group_idx <- picked[[1]]
      state_idx <- picked[[2]]
      chosen <- sample_from_group_state(group_idx, state_idx)
      if (!is.na(chosen)) {
        move_opinion(chosen, state_idx + 1L)
      }
    } else if (identical(event_type, "join")) {
      group_idx <- sample_weighted_state(join_term)
      if (totals[[group_idx]] < n) {
        chosen <- sample_outsider(group_idx, members, groups_of_individual, n)
        if (!is.na(chosen)) {
          add_membership(chosen, group_idx)
        }
      }
    } else if (identical(event_type, "leave_red")) {
      group_idx <- sample_weighted_state(leave_red_rate)
      state_pool <- camp_state_indices[[1]]
      state_weights <- state_counts[group_idx, state_pool]
      if (sum(state_weights) > 0L) {
        state_idx <- state_pool[[sample.int(length(state_pool), 1L, prob = state_weights)]]
        chosen <- sample_from_group_state(group_idx, state_idx)
        if (!is.na(chosen)) {
          remove_membership(chosen, group_idx)
        }
      }
    } else if (identical(event_type, "leave_blue")) {
      group_idx <- sample_weighted_state(leave_blue_rate)
      state_pool <- camp_state_indices[[2]]
      state_weights <- state_counts[group_idx, state_pool]
      if (sum(state_weights) > 0L) {
        state_idx <- state_pool[[sample.int(length(state_pool), 1L, prob = state_weights)]]
        chosen <- sample_from_group_state(group_idx, state_idx)
        if (!is.na(chosen)) {
          remove_membership(chosen, group_idx)
        }
      }
    } else if (identical(event_type, "infect_red")) {
      state_pool <- camp_state_indices[[1]]
      state_weights <- sir_counts[state_pool, 1L]
      if (sum(state_weights) > 0L) {
        state_idx <- state_pool[[sample.int(length(state_pool), 1L, prob = state_weights)]]
        chosen <- sample_from_compartment_state(1L, state_idx)
        if (!is.na(chosen)) {
          if (isTRUE(set_epidemic_state(chosen, 2L))) {
            record_infection(chosen, t)
          }
        }
      }
    } else if (identical(event_type, "infect_blue")) {
      state_pool <- camp_state_indices[[2]]
      state_weights <- sir_counts[state_pool, 1L]
      if (sum(state_weights) > 0L) {
        state_idx <- state_pool[[sample.int(length(state_pool), 1L, prob = state_weights)]]
        chosen <- sample_from_compartment_state(1L, state_idx)
        if (!is.na(chosen)) {
          if (isTRUE(set_epidemic_state(chosen, 2L))) {
            record_infection(chosen, t)
          }
        }
      }
    } else if (identical(event_type, "recover")) {
      if (sum(recovery_by_state) > 0L) {
        state_idx <- sample.int(length(recovery_by_state), 1L, prob = recovery_by_state)
        chosen <- sample_from_compartment_state(2L, state_idx)
        if (!is.na(chosen)) {
          set_epidemic_state(chosen, 3L)
        }
      }
    }

    event_counter <- event_counter + 1L
    if (event_counter %% record_every == 0L) {
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
    model_notes = "Exact individual/network simulation."
  )
}
