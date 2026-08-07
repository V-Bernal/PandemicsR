sample_social_event <- function(state, rates, params) {

  if (rates$social_rate <= 0)
    return(NULL)

  group_i <-
    if (params$m > 1)
      sample.int(
        params$m,
        1,
        prob = rates$lambda_i / rates$social_rate
      )
  else
    1

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

  if(length(state$outsiders[[group_i]] ) == 0){
    return(state)
  }

  outsiders <- state$outsiders
  members <- state$members
  opinions <- state$opinions

  # Join external to group
  chosen_idx <- sample.int(length(outsiders[[group_i]]),1)
  chosen <- outsiders[[group_i]][chosen_idx]
  outsiders[[group_i]][chosen_idx] <- outsiders[[group_i]][length(outsiders[[group_i]])]
  outsiders[[group_i]] <- outsiders[[group_i]][-length(outsiders[[group_i]])]

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

  state$outsiders <- outsiders
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
  idx_m <- which(state$members[[group_i]] == chosen)

  state$members[[group_i]][idx_m] <-
    state$members[[group_i]][length(state$members[[group_i]])]

  state$members[[group_i]] <-
    state$members[[group_i]][-length(state$members[[group_i]])]


  # Update individual membership
  state$groups_of_individual[[chosen]] <-
    setdiff(
      state$groups_of_individual[[chosen]],
      group_i
    )


  # Add to outsiders
  state$outsiders[[group_i]] <-
    c(state$outsiders[[group_i]], chosen)


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
  idx_m <- which(state$members[[group_i]] == chosen)

  state$members[[group_i]][idx_m] <-
    state$members[[group_i]][length(state$members[[group_i]])]

  state$members[[group_i]] <-
    state$members[[group_i]][-length(state$members[[group_i]])]


  # Update individual membership
  state$groups_of_individual[[chosen]] <-
    setdiff(
      state$groups_of_individual[[chosen]],
      group_i
    )


  # Add to outsiders
  state$outsiders[[group_i]] <-
    c(state$outsiders[[group_i]], chosen)


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

  members_g <- state$members[[group_i]]

  extreme_red <-
    members_g[state$opinions[members_g] == -2]


  if(length(extreme_red) == 0)
    return(state)


  chosen <- sample(extreme_red, 1)


  # Remove from members
  idx_m <- which(state$members[[group_i]] == chosen)

  state$members[[group_i]][idx_m] <-
    state$members[[group_i]][length(state$members[[group_i]])]

  state$members[[group_i]] <-
    state$members[[group_i]][-length(state$members[[group_i]])]


  # Update individual membership
  state$groups_of_individual[[chosen]] <-
    setdiff(
      state$groups_of_individual[[chosen]],
      group_i
    )


  # Add to outsiders
  state$outsiders[[group_i]] <-
    c(state$outsiders[[group_i]], chosen)


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

  if(length(state$Ei_blue[group_i]) == 0)
    return(state)

  members_g <- state$members[[group_i]]

  extreme_blue <-
    members_g[state$opinions[members_g] == 2]

  chosen <- sample(extreme_blue, 1)


  # Remove from members
  idx_m <- which(state$members[[group_i]] == chosen)

  state$members[[group_i]][idx_m] <-
    state$members[[group_i]][length(state$members[[group_i]])]

  state$members[[group_i]] <-
    state$members[[group_i]][-length(state$members[[group_i]])]


  # Update individual membership
  state$groups_of_individual[[chosen]] <-
    setdiff(
      state$groups_of_individual[[chosen]],
      group_i
    )


  # Add to outsiders
  state$outsiders[[group_i]] <-
    c(state$outsiders[[group_i]], chosen)


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
red_to_blue <- function(state, group_i){

  if(state$Ri[group_i] == 0 || state$Bi[group_i] == 0)
    return(state)

  chosen <- sample(state$red_members[[group_i]], 1)

  state$opinions[chosen] <- 1


  # Update all groups of individual
  for(g in state$groups_of_individual[[chosen]]){

    idx <- which(state$red_members[[g]] == chosen)

    if(length(idx) > 0){

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

  if (!is.null(state$epi) && state$epi[chosen] == state$S) {

    # Remove from susceptible red
    state$S_red_nodes <-
      state$S_red_nodes[state$S_red_nodes != chosen]

    # Add to susceptible blue
    state$S_blue_nodes <-
      c(state$S_blue_nodes, chosen)
  }

  # tracker warm up
  state$opinion_changes <-
    state$opinion_changes + 1

  return(state)
}



#---------------------------------------------------------
# Move 2: Moderate Blue -> Moderate Red
#---------------------------------------------------------
blue_to_red <- function(state, group_i){

  if(state$Ri[group_i] == 0 || state$Bi[group_i] == 0)
    return(state)

  chosen <- sample(state$blue_members[[group_i]], 1)

  state$opinions[chosen] <- -1


  for(g in state$groups_of_individual[[chosen]]){

    idx <- which(state$blue_members[[g]] == chosen)

    if(length(idx) > 0){

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

  if (!is.null(state$epi) && state$epi[chosen] == state$S) {

      # Remove from susceptible blue
      state$S_blue_nodes <-
        state$S_blue_nodes[state$S_blue_nodes != chosen]

      # Add to susceptible red
      state$S_red_nodes <-
        c(state$S_red_nodes, chosen)
    }

  # tracker warm up
  state$opinion_changes <-
    state$opinion_changes + 1

  return(state)
}



#---------------------------------------------------------
# Move 6: Moderate Red -> Extreme Red
#---------------------------------------------------------
radicalize_red <- function(state, group_i){

  if(state$Ri[group_i] == 0 ||
     state$Ei_red[group_i] == 0)
    return(state)

  candidates <-
    state$members[[group_i]][
      state$opinions[state$members[[group_i]]] == -1
    ]

  if(length(candidates)==0)
    return(state)


  chosen <- sample(candidates,1)

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

  if(state$Bi[group_i] == 0 ||
     state$Ei_blue[group_i] == 0)
    return(state)

  candidates <-
    state$members[[group_i]][
      state$opinions[state$members[[group_i]]] == 1
    ]

  if(length(candidates)==0)
    return(state)


  chosen <- sample(candidates,1)

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

  if(state$Ei_red[group_i] == 0 || state$Ri[group_i] > 0)
    return(state)


  candidates <-
    state$members[[group_i]][
      state$opinions[state$members[[group_i]]] == -2
    ]

  if(length(candidates)==0)
    return(state)


  chosen <- sample(candidates,1)

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


  chosen <- sample(candidates,1)

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
