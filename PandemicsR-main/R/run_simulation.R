#' Run the simulation
#'
#' @param params Simulation parameters.
#' @export
run_simulation <- function(params) {

  #==========================
  # 1. Initialization
  #==========================

  # Network
  network_state <- init_network(params)

  # Opinions
  opinion_state <- init_opinions(params)

  # Groups
  group_state <- init_groups(
    network_state$B0,
    opinion_state$opinions,
    params
  )

  #==========================
  # Build global state
  #==========================
  state <- list(

    opinions = opinion_state$opinions,

    members = group_state$members,
    red_members = group_state$red_members,
    blue_members = group_state$blue_members,
    outsiders = group_state$outsiders,
    groups_of_individual = group_state$groups_of_individual,

    Ri = group_state$Ri,
    Bi = group_state$Bi,

    Ei_red = group_state$Ei_red,
    Ei_blue = group_state$Ei_blue,


    B = network_state$B0,
    RIG = network_state$RIG0,
    rig_dirty = FALSE,

    members0 = group_state$members,
    in_group = logical(params$n)
  )

  #===========
  # Epidemic (optional)
  epi_state <- NULL

  if (isTRUE(params$runEpidemic)) {

    epi_state <- init_epidemic(
      params,
      opinion_state$opinions
    )

    state <- c(
      state,
      epi_state
    )

  }

  #===========
  # Stability if only epidemics start it immideitely
  #===========
  state$epidemic_started <- FALSE
  state$epidemic_started <- isTRUE(params$runEpidemic) && !isTRUE(params$runVoter) &&
      !isTRUE(params$runSchelling)


  stability_monitor <- create_stability_monitor(
    window = params$stability_window,
    threshold = params$stability_threshold
  )

  #===========
  # Trackers
  #===========
  trackers <- init_trackers(
    params,
    state$opinions
  )

  #==========================
  # 2. Gillespie loop
  #==========================

  # Gillespie algorithm
  start_time <- Sys.time()
  t <- 0
  stop_reason <- "Reached maximum time"

  while (t < params$t_max) {
    print(t)

    # Stopping criteria by full polarization, or ended epidemic.
    stop <- check_stopping(state, params)

    if (stop$stop) {
      stop_reason <- stop$reason
      break
    }

    # Gillespie step
    step <- gillespie_step(
      state,
      params,
      t
    )

    state <- step$state
    t <- step$t

    # stop no event (zero or negative rates), or maximum time
    if (step$stop) {
      stop_reason <- step$reason
      break
    }

    # Check stability

    # Start epidemic by time
    #if (params$epi_trigger == "time") {
    #  print('start epi time')
    #  stability_monitor$is_stable <- (t >= params$epi_time)
    #}else{
      stability_monitor <- update_stability(
        stability_monitor,
        state,
        time = t
      )
    #}

    if (
      !state$epidemic_started &&
      isTRUE(params$runEpidemic) &&
      isTRUE(stability_monitor$is_stable)
    ) {

      state$epidemic_started <- TRUE
      print('initialize epidemics')
      #update epidemic status
      epi_state <- init_epidemic(
        params,
        state$opinions
      )

      state[names(epi_state)] <- epi_state


        }

    #=======

    # ---- record step ----
    trackers$event_counter <- trackers$event_counter + 1

    if (t - trackers$last_record_time >= trackers$record_every) {

      trackers$last_record_time <- t

      trackers <- update_trackers(
        trackers,
        state,
        t,
        params
      )

      }
  }

  # end gillespie
  end_time <- Sys.time()
  comp_time <- end_time - start_time

  #==========================
  # Output
  #==========================
    result <- finalize_output(
      state,
      trackers,
      params,
      network_state,
      stop_reason,
      t,
      comp_time,
      trackers$event_counter
    )
  print(result )

    return(result)
}

