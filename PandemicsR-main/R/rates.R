compute_rates <- function(state, params){

  enabled_moves <- integer(0)

  if (params$runVoter)
    enabled_moves <- c(enabled_moves, 1, 2)

  if (params$runSchelling)
    enabled_moves <- c(enabled_moves, 3, 4, 5)

  if (params$runVoter && params$num_opinions == 4)
    enabled_moves <- c(enabled_moves, 6, 7, 8, 9)

  if (params$runSchelling && params$num_opinions == 4)
    enabled_moves <- c(enabled_moves, 10, 11)

  #==========================
  # Voter's rate of interaction among moderates
  #==========================
  voter_term <- numeric(params$m)

  Tot <- state$Ri + state$Bi# moderate_total
  Tot_g <- sapply(state$members, length) # group_size

  if (isTRUE(params$runVoter)) {
    valid <- (state$Ri > 0 & state$Bi > 0) #Tot >= 2
    voter_term[valid] <- params$gamma * (state$Ri[valid] * state$Bi[valid] / Tot[valid])
  }

  #==========================
  # Schelling rates
  ##==========================
  join_term <- numeric(params$m)
  leaveR_rate <- numeric(params$m)
  leaveB_rate <- numeric(params$m)

  if (isTRUE(params$runSchelling)) {
    # joining a group
    join_term <- params$c_param * sapply(state$outsiders, length) * sapply(state$members, length)

    if (isTRUE(params$scaled_n)) {
      join_term <- join_term/params$n
    }
    if (isTRUE(params$scaled_m)) {
      join_term <- join_term/params$m
    }

    # leaving a group rates
    # The rate of leaving a group is B+ or B- depending on the threshold
    frac_red <- ifelse(Tot > 0, state$Ri / Tot, 0)
    frac_blue <- ifelse(Tot > 0, state$Bi / Tot, 0)

    leaveR_rate <- ifelse(frac_red < params$T_threshold, params$beta_plus * state$Ri, params$beta_minus * state$Ri)
    leaveB_rate <- ifelse(frac_blue < params$T_threshold, params$beta_plus * state$Bi, params$beta_minus * state$Bi)

  }

  #==========================
  # New: de-radicalize rates
  #==========================
  rate_radicalize_red  <- numeric(params$m)
  rate_radicalize_blue <- numeric(params$m)
  rate_deradicalize_red <- numeric(params$m)
  rate_deradicalize_blue <- numeric(params$m)

  if (params$runVoter && params$num_opinions == 4) {

    factor <- ifelse(Tot_g > 0, 1 / Tot_g, 0)

    rate_radicalize_red  <- (params$alpha0 + params$alpha * state$Ei_red * factor ) * state$Ri
    rate_radicalize_blue <- (params$alpha0 + params$alpha * state$Ei_blue * factor )* state$Bi

    rate_deradicalize_red  <- (params$alpha0 + params$alpha_deradicalization * state$Ei_red * factor) * state$Ri
    rate_deradicalize_blue <- (params$alpha0 + params$alpha_deradicalization * state$Ei_blue * factor) * state$Bi

  }

  #==========================
  # New: leave rates radicals
  #==========================

  leaveR_rate_extreme <- numeric(params$m)
  leaveB_rate_extreme <- numeric(params$m)

  if (params$num_opinions == 4 && isTRUE(params$runSchelling)) {

    frac_red_e  <- ifelse(Tot_g > 0, state$Ei_red  / Tot_g, 0)
    frac_blue_e <- ifelse(Tot_g > 0, state$Ei_blue / Tot_g, 0)

    leaveR_rate_extreme <-
      ifelse(frac_red_e < params$T_threshold,
             params$beta_plus * state$Ei_red,
             params$beta_minus * state$Ei_red)

    leaveB_rate_extreme <-
      ifelse(frac_blue_e < params$T_threshold,
             params$beta_plus * state$Ei_blue,
             params$beta_minus * state$Ei_blue)

  }

  #==========================
  # Epidemic rates (NEW)
  #==========================
  infection_rate_tot <- 0
  recovery_rate_tot <- 0
  #beta_vec <- NULL
  epi_rate <- 0
  infection_red_rate <- 0
  infection_blue_rate <- 0

  if (isTRUE(params$runEpidemic)) {

    #I_count <- sum(state$epi == I)
    #prevalence <- state$I_count / params$n # global = random mixing

    # Opinion camp-based infection rates
    #beta_vec <- ifelse(state$opinions < 0, params$beta_red, params$beta_blue)

    # Aggregate epidemic rates

    # Infection: only for S individuals global random mixing.
    #infection_rate_tot <- (state$I_count / params$n) * sum(beta_vec[state$S_nodes])
    infection_red_rate <-
      params$beta_red *
      length(state$S_red_nodes) *
      state$I_count / params$n

    infection_blue_rate <-
      params$beta_blue *
      length(state$S_blue_nodes) *
      state$I_count / params$n

    infection_rate_tot <-
      infection_red_rate +
      infection_blue_rate


    # Recovery: only for I individuals
    recovery_rate_tot  <- params$gamma_epi * state$I_count

    epi_rate <- infection_rate_tot + recovery_rate_tot

  }

  #==========================
  # Global state rate
  #==========================
  # Each group has a params$lambda state rate
  lambda_i <- voter_term + join_term + leaveR_rate + leaveB_rate +
    rate_radicalize_red + rate_radicalize_blue +
    rate_deradicalize_red + rate_deradicalize_blue +
    leaveR_rate_extreme + leaveB_rate_extreme

  social_rate <- sum(lambda_i)

  lambda_tot <- social_rate + epi_rate

  list(
    voter_term = voter_term,
    join_term = join_term,
    leaveR_rate = leaveR_rate,
    leaveB_rate = leaveB_rate,

    rate_radicalize_red = rate_radicalize_red,
    rate_radicalize_blue = rate_radicalize_blue,
    rate_deradicalize_red = rate_deradicalize_red,
    rate_deradicalize_blue = rate_deradicalize_blue,

    leaveR_rate_extreme = leaveR_rate_extreme,
    leaveB_rate_extreme = leaveB_rate_extreme,

    infection_red_rate = infection_red_rate,
    infection_blue_rate = infection_blue_rate,

    infection_rate_tot = infection_rate_tot,
    recovery_rate_tot = recovery_rate_tot,

    lambda_i = lambda_i,
    social_rate = social_rate,
    epi_rate = epi_rate,
    lambda_tot = lambda_tot,

    enabled_moves = enabled_moves#,
    #beta_vec = beta_vec
  )
}

