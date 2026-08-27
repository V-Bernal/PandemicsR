#' Gillespie step
#'
#' @param params Simulation parameters.
#' @state state
#' @t time
#' @export
gillespie_step <- function(state, params, t){

  #==========================
  # Compute rates
  #==========================
  rates <- compute_rates(state, params)

  #print(rates)

  if (isTRUE(rates$lambda_tot <= 0)) {
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

  stopifnot(
    all(is.finite(rates$lambda_i)),
    all(rates$lambda_i >= 0),
    is.finite(rates$social_rate),
    rates$social_rate >= 0
  )

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
    if( isTRUE(state$epidemic_started) && isTRUE(params$runEpidemic) ){

      #if (!is.null(state$epi)) {

      #   stopifnot(
      #     length(state$epi) == length(state$opinions),
      #     length(state$epi) == length(state$groups_of_individual)
      #   )
      # }

      # cat("BEFORE epidemic event\n")
      # stopifnot(!anyNA(state$opinions))

      state <- apply_epidemic_event(
        state,
        rates,
        params,
        t
      )

      #cat("AFTER epidemic event\n")
      #stopifnot(!anyNA(state$opinions))

    }

  }

  return(list(
    state = state,
    t = t,
    stop = FALSE
  ))

}
