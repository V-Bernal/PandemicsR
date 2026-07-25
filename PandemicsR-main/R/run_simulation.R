#' Run the simulation
#'
#' @param params Simulation parameters.
#' @export
run_simulation <- function(params) {

  # 1. Initialization

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

  # Build global state
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

  # 2. Gillespie loop

    # real time
    start_time <- Sys.time()

    #===============
    # Progress bar
    #===============
    # withProgress(
    #   message = "Running simulation",
    #   detail = "Starting...",
    #   value = 0,
       #{

        # Gillespie algorithm
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

          # ---- record step ----
          trackers$event_counter <- trackers$event_counter + 1

          if (t - trackers$last_record_time >= trackers$record_every) {

            trackers$last_record_time <- t

            # setProgress(
            #   value = min(t / params$t_max, 1),
            #   detail = sprintf(
            #     "t = %.2f / %.2f | events = %d",
            #     t,
            #     params$t_max,
            #     trackers$event_counter
            #   )
            # )

            trackers <- update_trackers(
              trackers,
              state,
              t,
              params
            )

            }
        }
        # end gillespie
      #}
  #)
  end_time <- Sys.time()
  comp_time <- end_time - start_time

  # tracker_results <- finalize_trackers(
  #   trackers,
  #   params)


    # Output
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

