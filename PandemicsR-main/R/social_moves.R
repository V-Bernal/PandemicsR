#' sample_social_event
#'
#' @param params Simulation parameters.
#' @param state
#' @param rates
#' @export
sample_social_event <- function(state, rates, params) {

  if (rates$social_rate <= 0)
    return(NULL)

  group_i <-
    if (params$m > 1){
      sample.int(
        params$m,
        1,
        prob = rates$lambda_i / rates$social_rate)} else{1}

  Ri_g <- state$Ri[group_i]
  Bi_g <- state$Bi[group_i]

  voter_g  <- rates$voter_term[group_i]
  join_g   <- rates$join_term[group_i]
  leaveR_g <- rates$leaveR_rate[group_i]
  leaveB_g <- rates$leaveB_rate[group_i]

  rates_vec <- c(
    ifelse(voter_g > 0, voter_g/2, 0),
    ifelse(voter_g > 0, voter_g/2, 0),
    join_g,
    leaveR_g,
    leaveB_g,
    rates$rate_radicalize_red[group_i],
    rates$rate_radicalize_blue[group_i],
    rates$rate_deradicalize_red[group_i],
    rates$rate_deradicalize_blue[group_i],
    rates$leaveR_rate_extreme[group_i],
    rates$leaveB_rate_extreme[group_i]
  )

  rates_sub <- rates_vec[rates$enabled_moves]

  if (length(rates_sub) == 0 || sum(rates_sub) <= 0)
    return(NULL)

  move <- sample(
    rates$enabled_moves,
    1,
    prob = rates_sub / sum(rates_sub)
  )

  list(
    group = group_i,
    move  = move
  )
}

#---------------------------------------------------------
# apply_social_mov
#---------------------------------------------------------
apply_social_move <- function(state, move, group_i, params) {

  stopifnot(
    length(group_i) == 1,
    !is.na(group_i),
    group_i >= 1,
    group_i <= length(state$Ri),
    group_i <= length(state$Bi),
    group_i <= length(state$red_members),
    group_i <= length(state$blue_members)
  )

  if (move == 1) {

    state <- red_to_blue(state, group_i)

  } else if (move == 2) {

    state <- blue_to_red(state, group_i)

  } else if (move == 3) {

    state <- join_group(state, group_i, params)

  } else if (move == 4) {

    state <- leave_red(state, group_i)

  } else if (move == 5) {

    state <- leave_blue(state, group_i)

  } else if (move == 6) {

    state <- radicalize_red(state, group_i)

  } else if (move == 7) {

    state <- radicalize_blue(state, group_i)

  } else if (move == 8) {

    state <- deradicalize_red(state, group_i)

  } else if (move == 9) {

    state <- deradicalize_blue(state, group_i)

  } else if (move == 10) {

    state <- leave_extreme_red(state, group_i)

  } else if (move == 11) {

    state <- leave_extreme_blue(state, group_i)

  }

  return(state)
}


#---------------------------------------------------------
# move 3
#---------------------------------------------------------
join_group <- function(state, group_i, params){

  n_members <- length(state$members[[group_i]])
  if (n_members >= params$n) {
    return(state)
  }

  members <- state$members
  opinions <- state$opinions

  # Sample non-member node (rejection sampling or setdiff)
  if (n_members > params$n * 0.5) {
    outsiders <- setdiff(seq_len(params$n), members[[group_i]])
    if (length(outsiders) == 0) return(state)
    chosen <- outsiders[sample.int(length(outsiders), 1)]
  } else {
    repeat {
      chosen <- sample.int(params$n, 1)
      if (!(chosen %in% members[[group_i]])) break
    }
  }

  members[[group_i]] <- c(members[[group_i]], chosen)
  state$groups_of_individual[[chosen]] <- c(state$groups_of_individual[[chosen]], group_i)

  if (opinions[chosen] == -2) {
    state$Ei_red[group_i] <- state$Ei_red[group_i] + 1
  }

  if (opinions[chosen] == 2) {
    state$Ei_blue[group_i] <- state$Ei_blue[group_i] + 1
  }

  if (opinions[chosen] == -1) {
    state$red_members[[group_i]] <- c(state$red_members[[group_i]], chosen)
    state$Ri[group_i] <- state$Ri[group_i] + 1
  }

  if (opinions[chosen] == 1){
    state$blue_members[[group_i]] <- c(state$blue_members[[group_i]], chosen)
    state$Bi[group_i] <- state$Bi[group_i] + 1
  }

  state$members <- members
  state$rig_dirty <- TRUE

  state$in_group[chosen] <- TRUE

  # tracker warm up
  state$membership_changes <-
    state$membership_changes + 1

  return(state)
}

#---------------------------------------------------------
# move 4 Remove internal moderate Red
#---------------------------------------------------------
leave_red<- function(state, group_i, params){

  red_members_g <- state$red_members[[group_i]]

  if(length(red_members_g) == 0)
    return(state)

  chosen_idx <- sample.int(length(red_members_g), 1)
  chosen <- red_members_g[chosen_idx]


  # Remove from red_members
  state$red_members[[group_i]][chosen_idx] <-
    state$red_members[[group_i]][length(state$red_members[[group_i]])]

  state$red_members[[group_i]] <-
    state$red_members[[group_i]][-length(state$red_members[[group_i]])]


  # Remove from members
  idx_m <- match(chosen, state$members[[group_i]])
  if (!is.na(idx_m)) {
    len_m <- length(state$members[[group_i]])
    state$members[[group_i]][idx_m] <- state$members[[group_i]][len_m]
    state$members[[group_i]] <- state$members[[group_i]][-len_m]
  }


  # Update individual membership
  state$groups_of_individual[[chosen]] <-
    setdiff(
      state$groups_of_individual[[chosen]],
      group_i
    )


  # Update counters
  state$Ri[group_i] <-
    state$Ri[group_i] - 1


  state$rig_dirty <- TRUE

  state$in_group[chosen] <-
    length(state$groups_of_individual[[chosen]]) > 0

  # tracker warm up
  state$membership_changes <-
    state$membership_changes + 1

  return(state)
}


#---------------------------------------------------------
# Move 5: Remove internal moderate Blue
#---------------------------------------------------------
leave_blue <- function(state, group_i){

  blue_members_g <- state$blue_members[[group_i]]

  if(length(blue_members_g) == 0)
    return(state)

  chosen_idx <- sample.int(length(blue_members_g), 1)
  chosen <- blue_members_g[chosen_idx]


  # Remove from blue_members
  state$blue_members[[group_i]][chosen_idx] <-
    state$blue_members[[group_i]][length(state$blue_members[[group_i]])]

  state$blue_members[[group_i]] <-
    state$blue_members[[group_i]][-length(state$blue_members[[group_i]])]


  # Remove from members
  idx_m <- match(chosen, state$members[[group_i]])
  if (!is.na(idx_m)) {
    len_m <- length(state$members[[group_i]])
    state$members[[group_i]][idx_m] <- state$members[[group_i]][len_m]
    state$members[[group_i]] <- state$members[[group_i]][-len_m]
  }


  # Update individual membership
  state$groups_of_individual[[chosen]] <-
    setdiff(
      state$groups_of_individual[[chosen]],
      group_i
    )


  # Update counters
  state$Bi[group_i] <-
    state$Bi[group_i] - 1


  state$rig_dirty <- TRUE

  state$in_group[chosen] <-
    length(state$groups_of_individual[[chosen]]) > 0

  # tracker warm up
  state$membership_changes <-
    state$membership_changes + 1

  return(state)
}



#---------------------------------------------------------
# Move 10: Remove internal extreme Red
#---------------------------------------------------------
leave_extreme_red <- function(state, group_i){

  if(state$Ei_red[group_i] == 0)
    return(state)

  members_g <- state$members[[group_i]]

  extreme_red <-
    members_g[state$opinions[members_g] == -2]

  chosen_idx <- sample.int(length(extreme_red), 1)
  chosen <- extreme_red[chosen_idx]

  # Remove from members
  idx_m <- match(chosen, state$members[[group_i]])
  if (!is.na(idx_m)) {
    len_m <- length(state$members[[group_i]])
    state$members[[group_i]][idx_m] <- state$members[[group_i]][len_m]
    state$members[[group_i]] <- state$members[[group_i]][-len_m]
  }


  # Update individual membership
  state$groups_of_individual[[chosen]] <-
    setdiff(
      state$groups_of_individual[[chosen]],
      group_i
    )


  # Update extreme counter
  state$Ei_red[group_i] <-
    state$Ei_red[group_i] - 1


  state$rig_dirty <- TRUE

  state$in_group[chosen] <-
    length(state$groups_of_individual[[chosen]]) > 0

  # tracker warm up
  state$membership_changes <-
    state$membership_changes + 1

  return(state)
}



#---------------------------------------------------------
# Move 11: Remove internal extreme Blue
#---------------------------------------------------------
leave_extreme_blue <- function(state, group_i){

  if(state$Ei_blue[group_i] == 0)
    return(state)

  members_g <- state$members[[group_i]]

  extreme_blue <-
    members_g[state$opinions[members_g] == 2]

  chosen_idx <- sample.int(length(extreme_blue), 1)
  chosen <- extreme_blue[chosen_idx]


  # Remove from members
  idx_m <- match(chosen, state$members[[group_i]])
  if (!is.na(idx_m)) {
    len_m <- length(state$members[[group_i]])
    state$members[[group_i]][idx_m] <- state$members[[group_i]][len_m]
    state$members[[group_i]] <- state$members[[group_i]][-len_m]
  }


  # Update individual membership
  state$groups_of_individual[[chosen]] <-
    setdiff(
      state$groups_of_individual[[chosen]],
      group_i
    )


  # Update extreme counter
  state$Ei_blue[group_i] <-
    state$Ei_blue[group_i] - 1


  state$rig_dirty <- TRUE

  state$in_group[chosen] <-
    length(state$groups_of_individual[[chosen]]) > 0

  # tracker warm up
  state$membership_changes <-
    state$membership_changes + 1

  return(state)
}


#=========================================================
# Voter and radicalization events
#=========================================================

#---------------------------------------------------------
# Move 1: Moderate Red -> Moderate Blue
#---------------------------------------------------------
red_to_blue <- function(state, group_i) {

  if (state$Ri[group_i] == 0 || state$Bi[group_i] == 0)
    return(state)

  chosen_idx <- sample.int(length(state$red_members[[group_i]]), 1)
  chosen <- state$red_members[[group_i]][chosen_idx]
  #chosen <- sample(state$red_members[[group_i]], 1)

  stopifnot(
    length(chosen) == 1,
    !is.na(chosen),
    chosen >= 1,
    chosen <= length(state$opinions)
  )

  # Change opinion
  state$opinions[chosen] <- 1

  #-------------------------------------------------------
  # Update all groups of individual
  #-------------------------------------------------------
  for (g in state$groups_of_individual[[chosen]]) {

    idx <- which(state$red_members[[g]] == chosen)

    if (length(idx) > 0) {

      state$red_members[[g]][idx] <-
        state$red_members[[g]][length(state$red_members[[g]])]

      state$red_members[[g]] <-
        state$red_members[[g]][-length(state$red_members[[g]])]

      state$blue_members[[g]] <-
        c(state$blue_members[[g]], chosen)

      state$Ri[g] <- state$Ri[g] - 1
      state$Bi[g] <- state$Bi[g] + 1
    }
  }

  #-------------------------------------------------------
  # Update epidemic camp
  #-------------------------------------------------------
  if (isTRUE(state$epidemic_started)) {

    stopifnot(
      !is.null(state$epi),
      length(state$epi) == length(state$opinions),
      chosen >= 1,
      chosen <= length(state$epi)
    )

    if (state$epi[chosen] == state$S) {

      idx_sred <- match(chosen, state$S_red_nodes)
      if (!is.na(idx_sred)) {
        len_sred <- length(state$S_red_nodes)
        state$S_red_nodes[idx_sred] <- state$S_red_nodes[len_sred]
        state$S_red_nodes <- state$S_red_nodes[-len_sred]
      }

      state$S_blue_nodes <-
        c(state$S_blue_nodes, chosen)

      state$S_red <- state$S_red - 1
      state$S_blue <- state$S_blue + 1

    } else if (state$epi[chosen] == state$I) {

      idx_ired <- match(chosen, state$I_red_nodes)
      if (!is.na(idx_ired)) {
        len_ired <- length(state$I_red_nodes)
        state$I_red_nodes[idx_ired] <- state$I_red_nodes[len_ired]
        state$I_red_nodes <- state$I_red_nodes[-len_ired]
      }

      state$I_blue_nodes <-
        c(state$I_blue_nodes, chosen)

      state$I_red <- state$I_red - 1
      state$I_blue <- state$I_blue + 1

    } else if (state$epi[chosen] == state$R) {

      idx_rred <- match(chosen, state$R_red_nodes)
      if (!is.na(idx_rred)) {
        len_rred <- length(state$R_red_nodes)
        state$R_red_nodes[idx_rred] <- state$R_red_nodes[len_rred]
        state$R_red_nodes <- state$R_red_nodes[-len_rred]
      }

      state$R_blue_nodes <-
        c(state$R_blue_nodes, chosen)

      state$R_red <- state$R_red - 1
      state$R_blue <- state$R_blue + 1
    }

    state$total_red <- state$total_red - 1
    state$total_blue <- state$total_blue + 1
  }

  #-------------------------------------------------------
  # Tracker
  #-------------------------------------------------------
  state$opinion_changes <-
    state$opinion_changes + 1

  return(state)
}

#---------------------------------------------------------
# Move 2: Moderate Blue -> Moderate Red
#---------------------------------------------------------
blue_to_red <- function(state, group_i) {

  if (state$Ri[group_i] == 0 || state$Bi[group_i] == 0)
    return(state)

  chosen_idx <- sample.int(length(state$blue_members[[group_i]]), 1)
  chosen <- state$blue_members[[group_i]][chosen_idx]
  #chosen <- sample(state$blue_members[[group_i]], 1)

  # Change opinion
  state$opinions[chosen] <- -1

  #-------------------------------------------------------
  # Update all groups of individual
  #-------------------------------------------------------
  for (g in state$groups_of_individual[[chosen]]) {

    idx <- which(state$blue_members[[g]] == chosen)

    if (length(idx) > 0) {

      state$blue_members[[g]][idx] <-
        state$blue_members[[g]][length(state$blue_members[[g]])]

      state$blue_members[[g]] <-
        state$blue_members[[g]][-length(state$blue_members[[g]])]

      state$red_members[[g]] <-
        c(state$red_members[[g]], chosen)

      state$Bi[g] <- state$Bi[g] - 1
      state$Ri[g] <- state$Ri[g] + 1
    }
  }

  #-------------------------------------------------------
  # Update epidemic camp
  #-------------------------------------------------------
  if (isTRUE(state$epidemic_started)) {

    stopifnot(
      !is.null(state$epi),
      length(state$epi) == length(state$opinions),
      chosen >= 1,
      chosen <= length(state$epi)
    )

    if (state$epi[chosen] == state$S) {

      idx_sblue <- match(chosen, state$S_blue_nodes)
      if (!is.na(idx_sblue)) {
        len_sblue <- length(state$S_blue_nodes)
        state$S_blue_nodes[idx_sblue] <- state$S_blue_nodes[len_sblue]
        state$S_blue_nodes <- state$S_blue_nodes[-len_sblue]
      }

      state$S_red_nodes <-
        c(state$S_red_nodes, chosen)

      state$S_red <- state$S_red + 1
      state$S_blue <- state$S_blue - 1

    } else if (state$epi[chosen] == state$I) {

      idx_iblue <- match(chosen, state$I_blue_nodes)
      if (!is.na(idx_iblue)) {
        len_iblue <- length(state$I_blue_nodes)
        state$I_blue_nodes[idx_iblue] <- state$I_blue_nodes[len_iblue]
        state$I_blue_nodes <- state$I_blue_nodes[-len_iblue]
      }

      state$I_red_nodes <-
        c(state$I_red_nodes, chosen)

      state$I_red <- state$I_red + 1
      state$I_blue <- state$I_blue - 1

    } else if (state$epi[chosen] == state$R) {

      idx_rblue <- match(chosen, state$R_blue_nodes)
      if (!is.na(idx_rblue)) {
        len_rblue <- length(state$R_blue_nodes)
        state$R_blue_nodes[idx_rblue] <- state$R_blue_nodes[len_rblue]
        state$R_blue_nodes <- state$R_blue_nodes[-len_rblue]
      }

      state$R_red_nodes <-
        c(state$R_red_nodes, chosen)

      state$R_red <- state$R_red + 1
      state$R_blue <- state$R_blue - 1
    }

    state$total_red <- state$total_red - 1
    state$total_blue <- state$total_blue + 1
  }
  #-------------------------------------------------------
  # Tracker
  #-------------------------------------------------------
  state$opinion_changes <-
    state$opinion_changes + 1

  return(state)
}

#---------------------------------------------------------
# Move 6: Moderate Red -> Extreme Red
#---------------------------------------------------------
radicalize_red <- function(state, group_i){

  if(state$Ri[group_i] == 0)
    return(state)

  candidates <-
    state$members[[group_i]][
      state$opinions[state$members[[group_i]]] == -1
    ]

  if(length(candidates)==0)
    return(state)

  chosen_idx <- sample.int(length(candidates), 1)
  chosen <- candidates[chosen_idx]
  # chosen <- sample(candidates,1)

  state$opinions[chosen] <- -2


  for(g in state$groups_of_individual[[chosen]]){

    idx <- which(state$red_members[[g]] == chosen)

    if(length(idx)>0){

      state$red_members[[g]][idx] <-
        state$red_members[[g]][length(state$red_members[[g]])]

      state$red_members[[g]] <-
        state$red_members[[g]][-length(state$red_members[[g]])]

      state$Ri[g] <- state$Ri[g]-1
      state$Ei_red[g] <- state$Ei_red[g]+1
    }
  }

  # tracker warm up
  state$opinion_changes <-
    state$opinion_changes + 1

  return(state)
}



#---------------------------------------------------------
# Move 7: Moderate Blue -> Extreme Blue
#---------------------------------------------------------
radicalize_blue <- function(state, group_i){

  if(state$Bi[group_i] == 0 )
    return(state)

  candidates <-
    state$members[[group_i]][
      state$opinions[state$members[[group_i]]] == 1
    ]

  if(length(candidates)==0)
    return(state)

  chosen_idx <- sample.int(length(candidates), 1)
  chosen <- candidates[chosen_idx]
  #chosen <- sample(candidates,1)

  state$opinions[chosen] <- 2


  for(g in state$groups_of_individual[[chosen]]){

    idx <- which(state$blue_members[[g]] == chosen)

    if(length(idx)>0){

      state$blue_members[[g]][idx] <-
        state$blue_members[[g]][length(state$blue_members[[g]])]

      state$blue_members[[g]] <-
        state$blue_members[[g]][-length(state$blue_members[[g]])]

      state$Bi[g] <- state$Bi[g]-1
      state$Ei_blue[g] <- state$Ei_blue[g]+1
    }
  }

  # tracker warm up
  state$opinion_changes <-
    state$opinion_changes + 1

  return(state)
}



#---------------------------------------------------------
# Move 8: Extreme Red -> Moderate Red
#---------------------------------------------------------
deradicalize_red <- function(state, group_i){

  if(state$Ei_red[group_i] == 0 || state$Ri[group_i] == 0)
    return(state)


  candidates <-
    state$members[[group_i]][
      state$opinions[state$members[[group_i]]] == -2
    ]

  if(length(candidates)==0)
    return(state)

  chosen_idx <- sample.int(length(candidates), 1)
  chosen <- candidates[chosen_idx]
  #chosen <- sample(candidates,1)

  state$opinions[chosen] <- -1


  for(g in state$groups_of_individual[[chosen]]){

    state$red_members[[g]] <-
      c(state$red_members[[g]], chosen)

    state$Ri[g] <- state$Ri[g]+1
    state$Ei_red[g] <- state$Ei_red[g]-1
  }

  # tracker warm up
  state$opinion_changes <-
    state$opinion_changes + 1

  return(state)
}



#---------------------------------------------------------
# Move 9: Extreme Blue -> Moderate Blue
#---------------------------------------------------------
deradicalize_blue <- function(state, group_i){

  if(state$Ei_blue[group_i] == 0 || state$Bi[group_i] == 0)
    return(state)

  candidates <-
    state$members[[group_i]][
      state$opinions[state$members[[group_i]]] == 2
    ]

  if(length(candidates)==0)
    return(state)

  chosen_idx <- sample.int(length(candidates), 1)
  chosen <- candidates[chosen_idx]
  #chosen <- sample(candidates,1)

  state$opinions[chosen] <- 1


  for(g in state$groups_of_individual[[chosen]]){

    state$blue_members[[g]] <-
      c(state$blue_members[[g]], chosen)

    state$Bi[g] <- state$Bi[g]+1
    state$Ei_blue[g] <- state$Ei_blue[g]-1
  }

  # tracker warm up
  state$opinion_changes <-
    state$opinion_changes + 1

  return(state)
}

