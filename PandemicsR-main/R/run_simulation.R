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

  # Epidemic (optional)
  epi_state <- NULL

  if (isTRUE(params$runEpidemic)) {

    epi_state <- init_epidemic(
      params,
      opinion_state$opinions
    )

  }

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
  # add tracker of warm up
  # NEW: epidemic warm-up tracking
  state$events <- 0
  state$opinion_changes <- 0
  state$membership_changes <- 0
  state$stable <- FALSE
  state$stable_windows <- 0
  t0 <-0
  #===========

  # Add epidemic state if active
  if (!is.null(epi_state)) {

    state <- c(
      state,
      epi_state
    )

  }

  # Trackers
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

    # Stopping criteria
    #stop <- check_stopping(state, params)

    # if (stop$stop) {
    #   stop_reason <- stop$reason
    #   break
    # }

    # Gillespie step
    step <- gillespie_step(
      state,
      params,
      t
    )

    state <- step$state
    t <- step$t

    if(step$stop)
      break

    #======= check warm up track reset every record_every
    check_stable <- check_stability(state, t, t0)
    state <- check_stable[[1]]
    t0 <- check_stable[[2]]

    if(isTRUE(state$stable)){
      print('stability reached')
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

    return(result)
}

