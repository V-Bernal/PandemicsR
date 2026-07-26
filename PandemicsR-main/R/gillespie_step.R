gillespie_step <- function(state, params, t){

  #==========================
  # Compute rates
  #==========================
  rates <- compute_rates(state, params)

  if (rates$lambda_tot <= 0) {
    return(list(
      state = state,
      t = t,
      stop = TRUE,
      reason = "No events possible"
    ))
  }

  #==========================
  # Gillespie time
  #==========================
  dt <- stats::rexp(1, rates$lambda_tot)
  t <- t + dt

  if (t >= params$t_max) {
    return(list(
      state = state,
      t = t,
      stop = TRUE,
      reason = "Maximum time reached"
    ))
  }

  #==========================
  # Select Social or Epidemic event
  #==========================
  u_event <- stats::runif(1)

  if (u_event < rates$social_rate / rates$lambda_tot){

    event <- sample_social_event(
      state,
      rates,
      params
    )

    if (!is.null(event)){
      state <- apply_social_move(
        state,
        event$move,
        event$group,
        params
      )
    }

  } else {

    # Check epidemic activation
    if( t >= params$epi_time){
      state <- apply_epidemic_event(
        state,
        rates,
        params,
        t
      )
    }

  }

  return(list(
    state = state,
    t = t,
    stop = FALSE
  ))

}
