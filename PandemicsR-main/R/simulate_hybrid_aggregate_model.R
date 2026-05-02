#' Simulate a Fast Aggregate Voter-Schelling-SIR Model
#'
#' @author OpenAI Codex
#'
#' @name simulate_hybrid_aggregate_model
#'
#' @inheritParams simulate_hybrid_model
#' @param max_steps Maximum number of aggregate integration steps.
#'
#' @return A list with the same plotting summaries as \code{simulate_hybrid_model}.
#' @export
simulate_hybrid_aggregate_model <- function(
    n, m, t_max, lambda, c_param,
    gamma_light, gamma_dark = gamma_light,
    infected_dark_multiplier = 1,
    beta_plus, beta_minus, T_threshold,
    num_opinions = 2,
    run_voter = TRUE, run_schelling = TRUE,
    beta_red = 0.5, beta_blue = 0.25, gamma_sir = 0.2,
    initial_infected_fraction = 0.05,
    max_steps = 600L) {

  if (!(num_opinions %in% c(2L, 4L))) {
    stop("This implementation currently supports 2 or 4 opinion states.")
  }

  levels_vec <- get_levels_vec(num_opinions)
  state_labels <- get_state_labels(num_opinions)
  state_camps <- get_state_camp_index(num_opinions)
  dark_states <- get_dark_state_flags(num_opinions)
  n_states <- length(levels_vec)
  camp_labels <- get_camp_labels()
  camp_state_indices <- lapply(seq_along(camp_labels), function(idx) which(state_camps == idx))
  compartment_labels <- c("S", "I", "R")

  if (num_opinions <= 2L) {
    gamma_dark <- gamma_light
    infected_dark_multiplier <- 1
  }

  parse_initial_infected <- function(value, population_size) {
    if (!is.finite(value) || value <= 0) {
      return(0)
    }

    if (value > 1) {
      return(min(population_size, round(value)))
    }

    min(population_size, max(1, round(value * population_size)))
  }

  state_counts <- rep(floor(n / n_states), n_states)
  if (sum(state_counts) < n) {
    state_counts[seq_len(n - sum(state_counts))] <- state_counts[seq_len(n - sum(state_counts))] + 1
  }
  names(state_counts) <- state_labels

  initial_infected_n <- parse_initial_infected(initial_infected_fraction, n)
  initial_infected_by_state <- state_counts * initial_infected_n / n

  sir_counts <- matrix(
    0,
    nrow = n_states,
    ncol = length(compartment_labels),
    dimnames = list(state_labels, compartment_labels)
  )
  sir_counts[, "S"] <- state_counts - initial_infected_by_state
  sir_counts[, "I"] <- initial_infected_by_state

  base_state_rates <- rep(gamma_light, n_states)
  base_state_rates[dark_states] <- gamma_dark
  infected_state_rates <- base_state_rates
  infected_state_rates[dark_states] <- infected_state_rates[dark_states] * infected_dark_multiplier

  sum_state_fraction <- function(current_counts) {
    rowSums(current_counts) / n
  }

  sum_sir_counts <- function(current_counts) {
    c(
      S = sum(current_counts[, "S"]),
      I = sum(current_counts[, "I"]),
      R = sum(current_counts[, "R"])
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

  t_max <- max(0, as.numeric(t_max))
  n_steps <- if (t_max <= 0) 0L else min(max_steps, max(1L, ceiling(t_max / 0.25)))
  dt <- if (n_steps > 0L) t_max / n_steps else 0

  frac_history <- vector("list", n_steps + 1L)
  sir_history <- vector("list", n_steps + 1L)
  camp_sir_history <- vector("list", n_steps + 1L)
  time_history <- seq(0, t_max, length.out = n_steps + 1L)

  infection_events <- data.frame(
    time = rep(0, length(camp_labels)),
    camp = camp_labels,
    count = vapply(camp_state_indices, function(state_idx) {
      sum(initial_infected_by_state[state_idx])
    }, numeric(1)),
    stringsAsFactors = FALSE
  )

  record_snapshot <- function(idx) {
    frac_history[[idx]] <<- sum_state_fraction(sir_counts)
    sir_history[[idx]] <<- sum_sir_counts(sir_counts)
    camp_sir_history[[idx]] <<- sum_camp_sir_counts(sir_counts)
  }

  record_snapshot(1L)

  for (step_idx in seq_len(n_steps)) {
    current_time <- time_history[[step_idx + 1L]]

    if (isTRUE(run_voter)) {
      state_totals <- rowSums(sir_counts)
      cumulative <- cumsum(state_totals)
      lower_counts <- c(0, cumulative[-length(cumulative)])
      higher_counts <- n - cumulative
      delta <- matrix(0, nrow = n_states, ncol = length(compartment_labels))

      for (state_idx in seq_len(n_states)) {
        if (state_totals[[state_idx]] <= 0) {
          next
        }

        for (comp_idx in seq_along(compartment_labels)) {
          available <- sir_counts[state_idx, comp_idx]
          if (available <= 0) {
            next
          }

          changer_rate <- if (comp_idx == 2L) infected_state_rates[[state_idx]] else base_state_rates[[state_idx]]
          flow_down <- if (state_idx > 1L) dt * changer_rate * available * lower_counts[[state_idx]] / n else 0
          flow_up <- if (state_idx < n_states) dt * changer_rate * available * higher_counts[[state_idx]] / n else 0
          total_flow <- flow_down + flow_up

          if (total_flow > available && total_flow > 0) {
            scale <- available / total_flow
            flow_down <- flow_down * scale
            flow_up <- flow_up * scale
          }

          delta[state_idx, comp_idx] <- delta[state_idx, comp_idx] - flow_down - flow_up
          if (flow_down > 0) {
            delta[state_idx - 1L, comp_idx] <- delta[state_idx - 1L, comp_idx] + flow_down
          }
          if (flow_up > 0) {
            delta[state_idx + 1L, comp_idx] <- delta[state_idx + 1L, comp_idx] + flow_up
          }
        }
      }

      sir_counts <- sir_counts + delta
      sir_counts[sir_counts < 0] <- 0
    }

    infected_total <- sum(sir_counts[, "I"])
    new_infections_by_camp <- numeric(length(camp_labels))

    if (infected_total > 0) {
      for (camp_idx in seq_along(camp_labels)) {
        state_idx <- camp_state_indices[[camp_idx]]
        beta <- if (camp_idx == 1L) beta_red else beta_blue
        infection_rates <- beta * infected_total / n * sir_counts[state_idx, "S"]
        new_infections <- pmin(sir_counts[state_idx, "S"], dt * infection_rates)

        sir_counts[state_idx, "S"] <- sir_counts[state_idx, "S"] - new_infections
        sir_counts[state_idx, "I"] <- sir_counts[state_idx, "I"] + new_infections
        new_infections_by_camp[[camp_idx]] <- sum(new_infections)
      }
    }

    recoveries <- pmin(sir_counts[, "I"], dt * gamma_sir * sir_counts[, "I"])
    sir_counts[, "I"] <- sir_counts[, "I"] - recoveries
    sir_counts[, "R"] <- sir_counts[, "R"] + recoveries

    if (any(new_infections_by_camp > 0)) {
      infection_events <- rbind(
        infection_events,
        data.frame(
          time = rep(current_time, length(camp_labels)),
          camp = camp_labels,
          count = new_infections_by_camp,
          stringsAsFactors = FALSE
        )
      )
    }

    record_snapshot(step_idx + 1L)
  }

  frac_mat <- do.call(rbind, frac_history)
  colnames(frac_mat) <- state_labels

  sir_mat <- do.call(rbind, sir_history)
  colnames(sir_mat) <- compartment_labels

  final_camp_sir <- camp_sir_history[[length(camp_sir_history)]]
  final_camp_sizes <- rowSums(final_camp_sir)
  camp_attack_rate <- ifelse(final_camp_sizes > 0, 1 - final_camp_sir[, "S"] / final_camp_sizes, NA_real_)
  names(camp_attack_rate) <- camp_labels

  list(
    B0 = NULL,
    B = NULL,
    RIG0 = NULL,
    RIG = NULL,
    opinions = NULL,
    opinion_history = NULL,
    num_opinions = num_opinions,
    members0 = vector("list", m),
    members = vector("list", m),
    frac_mat = frac_mat,
    sir_mat = sir_mat,
    camp_sir_history = camp_sir_history,
    infection_events = infection_events[infection_events$count > 0, , drop = FALSE],
    time_history = time_history,
    sir_state = NULL,
    final_time = t_max,
    state_labels = state_labels,
    camp_labels = camp_labels,
    overall_attack_rate = 1 - sum(sir_counts[, "S"]) / n,
    camp_attack_rate = camp_attack_rate,
    final_camp_sir = final_camp_sir,
    simulation_mode = "aggregate",
    graph_available = FALSE
  )
}
