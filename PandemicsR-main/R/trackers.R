init_trackers <- function(params, opinions) {

  n <- params$n
  num_opinions <- params$num_opinions

  # Gillespie recording
  record_every <- params$t_max / 500
  max_records  <- ceiling(params$t_max / record_every) + 1

  # Opinion levels
  levels <- if (num_opinions == 4) {
    c(-2, -1, 1, 2)
  } else {
    c(-1, 1)
  }

  # Time history
  time_hist <- matrix(
    NA_real_,
    nrow = 1,
    ncol = max_records
  )
  time_hist[1] <- 0

  # Opinion fractions
  frac_mat <- matrix(
    NA_real_,
    nrow = length(levels),
    ncol = max_records
  )

  frac_mat[, 1] <- sapply(
    levels,
    function(op) sum(opinions == op) / n
  )

  opinion_history <- matrix(
    NA,
    nrow = params$n,
    ncol = 2
  )
  opinion_history[,1] <- opinions

  trackers <- list(

    # Event tracking
    event_counter    = 0,
    record_every     = record_every,
    max_records      = max_records,
    last_record_time = 0,
    record_idx       = 1,

    # Opinion tracking
    levels    = levels,
    time_hist = time_hist,
    frac_mat  = frac_mat,
    opinion_history = opinion_history

  )

  if (isTRUE(params$runEpidemic)) {

    trackers$S_hist <- rep(NA_real_, max_records)
    trackers$I_hist <- rep(NA_real_, max_records)
    trackers$R_hist <- rep(NA_real_, max_records)

    trackers$S_hist_red  <- rep(NA_real_, max_records)
    trackers$I_hist_red  <- rep(NA_real_, max_records)
    trackers$R_hist_red  <- rep(NA_real_, max_records)

    trackers$S_hist_blue <- rep(NA_real_, max_records)
    trackers$I_hist_blue <- rep(NA_real_, max_records)
    trackers$R_hist_blue <- rep(NA_real_, max_records)

    trackers$S_hist_grp <- rep(NA_real_, max_records)
    trackers$I_hist_grp <- rep(NA_real_, max_records)
    trackers$R_hist_grp <- rep(NA_real_, max_records)

    trackers$S_hist_out <- rep(NA_real_, max_records)
    trackers$I_hist_out <- rep(NA_real_, max_records)
    trackers$R_hist_out <- rep(NA_real_, max_records)
  }

  trackers
}

#-----------------------------
update_trackers <- function(trackers, state, t, params) {

  trackers$opinion_history[,2] <- state$opinions

  # Advance record index
  trackers$record_idx <- trackers$record_idx + 1
  idx <- trackers$record_idx

  #=========================
  # Time
  #=========================
  trackers$time_hist[idx] <- t

  #=========================
  # Opinion fractions
  #=========================
  trackers$frac_mat[, idx] <- sapply(
    trackers$levels,
    function(op) sum(state$opinions == op) / params$n
  )

  #=========================
  # Epidemic tracking
  #=========================
  if (isTRUE(params$runEpidemic)) {

    isS <- state$epi == state$S
    isI <- state$epi == state$I
    isR <- state$epi == state$R

    grp_idx <- state$in_group
    out_idx <- !grp_idx

    # Overall
    trackers$S_hist[idx] <- state$S_count
    trackers$I_hist[idx] <- state$I_count
    trackers$R_hist[idx] <- state$R_count

    # Opinion camps
    trackers$S_hist_red[idx]  <- state$S_red
    trackers$I_hist_red[idx]  <- state$I_red
    trackers$R_hist_red[idx]  <- state$R_red

    trackers$S_hist_blue[idx] <- state$S_blue
    trackers$I_hist_blue[idx] <- state$I_blue
    trackers$R_hist_blue[idx] <- state$R_blue

    # Membership
    trackers$S_hist_grp[idx] <- sum(isS & grp_idx)
    trackers$I_hist_grp[idx] <- sum(isI & grp_idx)
    trackers$R_hist_grp[idx] <- sum(isR & grp_idx)

    trackers$S_hist_out[idx] <- sum(isS & out_idx)
    trackers$I_hist_out[idx] <- sum(isI & out_idx)
    trackers$R_hist_out[idx] <- sum(isR & out_idx)
  }

  trackers
}

######################
finalize_trackers <- function(trackers, params) {

  #=========================
  # Finalize opinion history
  #=========================

  time_hist <- as.vector(t(trackers$time_hist))

  frac_mat <- t(trackers$frac_mat)

  colnames(frac_mat) <- trackers$levels

  # Remove unused records
  id_t <- which(!is.na(time_hist))

  time_hist <- time_hist[id_t]
  frac_mat <- frac_mat[id_t, , drop = FALSE]


  results <- list(
    time_hist = time_hist,
    frac_mat = frac_mat
  )


  #=========================
  # Epidemic histories
  #=========================

  if (isTRUE(params$runEpidemic)) {

    results$SIR_df <- data.frame(
      time = time_hist,
      S = trackers$S_hist[id_t] / params$n,
      I = trackers$I_hist[id_t] / params$n,
      R = trackers$R_hist[id_t] / params$n
    )


    results$SIR_df_opinion_red <- data.frame(
      time = time_hist,
      S = trackers$S_hist_red[id_t] / params$n,
      I = trackers$I_hist_red[id_t] / params$n,
      R = trackers$R_hist_red[id_t] / params$n
    )


    results$SIR_df_opinion_blue <- data.frame(
      time = time_hist,
      S = trackers$S_hist_blue[id_t] / params$n,
      I = trackers$I_hist_blue[id_t] / params$n,
      R = trackers$R_hist_blue[id_t] / params$n
    )


    results$SIR_df_grp <- data.frame(
      time = time_hist,
      S = trackers$S_hist_grp[id_t] / params$n,
      I = trackers$I_hist_grp[id_t] / params$n,
      R = trackers$R_hist_grp[id_t] / params$n
    )


    results$SIR_df_out <- data.frame(
      time = time_hist,
      S = trackers$S_hist_out[id_t] / params$n,
      I = trackers$I_hist_out[id_t] / params$n,
      R = trackers$R_hist_out[id_t] / params$n
    )

  } else {

    results$SIR_df <- NULL
    results$SIR_df_opinion_red <- NULL
    results$SIR_df_opinion_blue <- NULL
    results$SIR_df_grp <- NULL
    results$SIR_df_out <- NULL

  }


  results
}
