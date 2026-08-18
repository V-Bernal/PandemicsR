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

  # epidemic state
  # epi_state <- init_epidemic(
  #   params,
  #   opinion_state$opinions
  # )
  #
  # state <- c(
  #   state,
  #   epi_state
  # )

  state$epidemic_started <- FALSE
  stability_monitor <- create_stability_monitor(
    window = params$t_max/20,
    threshold = 0.01,
    persistence = 5
  )

  #===========

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

    if(step$stop)
      break

    # New Check stability
    stability_monitor <- update_stability(
      stability_monitor,
      state,
      time = t
    )

    # Start epidemic once social system is stable
    if (!state$epidemic_started && isTRUE(params$runEpidemic)){

      state$epidemic_started <- stability_monitor$is_stable

      if (state$epidemic_started) {

          epi_state <- init_epidemic(
            params,
            opinion_state$opinions
          )

          state <- c(
            state,
            epi_state
          )

        }
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

