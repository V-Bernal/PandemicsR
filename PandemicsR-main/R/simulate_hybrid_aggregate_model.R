#' Simulate a Fast Aggregate Voter-Schelling-SIR Model
#'
#' This is a large-population approximation using the same public parameters as
#' the main-compatible exact simulator. It preserves the dashboard API while
#' avoiding materializing individual networks for very large n.
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
    max_steps = 600L) {

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
  max_steps <- max(1L, as.integer(round(max_steps)))

  levels_vec <- get_levels_vec(num_opinions)
  state_labels <- get_state_labels(num_opinions)
  state_camps <- get_state_camp_index(num_opinions)
  n_states <- length(levels_vec)
  camp_labels <- get_camp_labels()
  camp_state_indices <- lapply(seq_along(camp_labels), function(idx) which(state_camps == idx))
  compartment_labels <- c("S", "I", "R")

  red_mod_idx <- match(-1, levels_vec)
  blue_mod_idx <- match(1, levels_vec)
  red_ext_idx <- match(-2, levels_vec)
  blue_ext_idx <- match(2, levels_vec)

  parse_initial_infected <- function(value, population_size) {
    if (!is.finite(value) || value <= 0) {
      return(0)
    }
    if (value > 1) {
      return(min(population_size, floor(value)))
    }
    min(population_size, floor(value * population_size))
  }

  state_counts <- rep(floor(n / n_states), n_states)
  if (sum(state_counts) < n) {
    state_counts[seq_len(n - sum(state_counts))] <- state_counts[seq_len(n - sum(state_counts))] + 1
  }
  names(state_counts) <- state_labels

  initial_infected_n <- if (isTRUE(run_epidemic)) parse_initial_infected(initial_infected_fraction, n) else 0
  initial_infected_by_state <- state_counts * initial_infected_n / n

  sir_counts <- matrix(
    0,
    nrow = n_states,
    ncol = length(compartment_labels),
    dimnames = list(state_labels, compartment_labels)
  )
  sir_counts[, "S"] <- state_counts - initial_infected_by_state
  sir_counts[, "I"] <- initial_infected_by_state

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

  normalize_population <- function() {
    total_count <- sum(sir_counts)
    if (is.finite(total_count) && total_count > 0 && abs(total_count - n) > 1e-7) {
      sir_counts <<- sir_counts * (n / total_count)
    }
  }

  record_snapshot(1L)

  for (step_idx in seq_len(n_steps)) {
    current_time <- time_history[[step_idx + 1L]]

    if (isTRUE(run_voter)) {
      delta <- matrix(0, nrow = n_states, ncol = length(compartment_labels))
      state_totals <- rowSums(sir_counts)

      if (!is.na(red_mod_idx) && !is.na(blue_mod_idx)) {
        red_total <- state_totals[[red_mod_idx]]
        blue_total <- state_totals[[blue_mod_idx]]
        moderate_total <- red_total + blue_total
        if (moderate_total > 0) {
          red_to_blue <- dt * gamma * red_total * blue_total / moderate_total / 2
          blue_to_red <- red_to_blue
          red_to_blue <- min(red_to_blue, red_total)
          blue_to_red <- min(blue_to_red, blue_total)
          if (red_total > 0 && red_to_blue > 0) {
            moved <- sir_counts[red_mod_idx, ] * (red_to_blue / red_total)
            delta[red_mod_idx, ] <- delta[red_mod_idx, ] - moved
            delta[blue_mod_idx, ] <- delta[blue_mod_idx, ] + moved
          }
          if (blue_total > 0 && blue_to_red > 0) {
            moved <- sir_counts[blue_mod_idx, ] * (blue_to_red / blue_total)
            delta[blue_mod_idx, ] <- delta[blue_mod_idx, ] - moved
            delta[red_mod_idx, ] <- delta[red_mod_idx, ] + moved
          }
        }
      }

      if (num_opinions == 4L) {
        state_totals <- rowSums(sir_counts + delta)
        total_group_proxy <- max(1, sum(state_totals))

        red_rad <- dt * (alpha0 + alpha * state_totals[[red_ext_idx]] / total_group_proxy) * state_totals[[red_mod_idx]]
        blue_rad <- dt * (alpha0 + alpha * state_totals[[blue_ext_idx]] / total_group_proxy) * state_totals[[blue_mod_idx]]
        red_derad <- dt * (alpha0 + alpha_deradicalization * state_totals[[red_ext_idx]] / total_group_proxy) * state_totals[[red_mod_idx]]
        blue_derad <- dt * (alpha0 + alpha_deradicalization * state_totals[[blue_ext_idx]] / total_group_proxy) * state_totals[[blue_mod_idx]]

        red_rad <- min(red_rad, state_totals[[red_mod_idx]])
        blue_rad <- min(blue_rad, state_totals[[blue_mod_idx]])
        red_derad <- min(red_derad, state_totals[[red_ext_idx]])
        blue_derad <- min(blue_derad, state_totals[[blue_ext_idx]])

        move_between_states <- function(from_idx, to_idx, amount) {
          if (amount <= 0 || state_totals[[from_idx]] <= 0) {
            return(invisible(FALSE))
          }
          moved <- (sir_counts[from_idx, ] + delta[from_idx, ]) * (amount / state_totals[[from_idx]])
          delta[from_idx, ] <<- delta[from_idx, ] - moved
          delta[to_idx, ] <<- delta[to_idx, ] + moved
          invisible(TRUE)
        }

        move_between_states(red_mod_idx, red_ext_idx, red_rad)
        move_between_states(blue_mod_idx, blue_ext_idx, blue_rad)
        move_between_states(red_ext_idx, red_mod_idx, red_derad)
        move_between_states(blue_ext_idx, blue_mod_idx, blue_derad)
      }

      sir_counts <- sir_counts + delta
      sir_counts[sir_counts < 0] <- 0
      normalize_population()
    }

    new_infections_by_camp <- numeric(length(camp_labels))
    if (isTRUE(run_epidemic)) {
      infected_total <- sum(sir_counts[, "I"])
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
      normalize_population()
    }

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
    graph_available = FALSE,
    population_size = n,
    group_count = m,
    event_count = NA_integer_,
    stop_reason = "Reached maximum time",
    model_notes = "Fast aggregate approximation using the same public parameters as main-compatible exact mode."
  )
}
