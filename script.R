library(igraph)
library(Matrix)
library(PandemicsR)
#==========================
# 1. User parameters
#==========================

n <- 15 # number of individuals
m <- 6 #
num_opinions <- 2

lambda <- 0.5 # "Init: Bipartite link density lambda"
c_param <- 0.81 # rate of joining a group
gamma <- 0.5 # Voter rate parameter

beta_plus <- 0.7 # leave rate
beta_minus <- 0.3 # leave rate
T_threshold <- 0.5 # if fraction of less than T

t_max <- 1000 # Gillespie time

#==========================
# 2. Model parameters
#==========================
# Homogeneous weights
ind_w <- rep(lambda, n) #ind_w <- runif(n, 1, 2)
grp_w <- rep(lambda * n / m, m) #grp_w <- runif(m, 1, 2) * n / m

# bipartite and RIG
B0 <- generate_bipartite(n, m, ind_w, grp_w) #B0 <- generate_bipartite(n, m, ind_w, grp_w, lambda)
RIG0 <- bipartite_to_rig(B0)
RIG <- RIG0; bipartite <- B0

# Voter tracker
opinions <- initialize_opinions_multi(n, num_opinions)
#opinions0 <- opinions
opinion_history <- matrix(opinions, ncol = 1)

levels <- sort(unique(as.vector(opinions)))
frac_mat <- sapply(levels, function(op) {
  sum(opinions == op)/n
})
  
# Group tracker
members <- vector("list", m)
red_members <- vector("list", m)
blue_members <- vector("list", m)
outsiders <- vector("list", m)
groups_of_individual <- vector("list", n)

for (g in seq_len(m)) {
  #group
  ids <- which(bipartite[, g] == 1)
  members[[g]] <- ids
  #opinion in group
  red_members[[g]] <- ids[opinions[ids] == -1]
  blue_members[[g]] <- ids[opinions[ids] == +1]
  outsiders[[g]] <- setdiff(seq_len(n), ids)
  for (i in ids) groups_of_individual[[i]] <- c(groups_of_individual[[i]], g)
}

Ri <- sapply(red_members, length)
Bi <- sapply(blue_members, length)
rig_dirty <- FALSE

#----- NEW
event_counter <- 0
record_every <- n
members0<-members

#==========================
# Combined Voter and Schelling model dynamics Gillespie algorithm
#==========================
t <- 0
while (t < t_max) {

  # The rate of interaction
  Tot <- Ri + Bi
  #voter_term <- gamma * (Ri * Bi/ Tot) # Contact rate ind vd ind: gamma*k
  #voter_term[Tot < 2] <- 0
  
  # new
  voter_term <- numeric(m)
  valid <- Tot >= 2
  voter_term[valid] <- gamma * (Ri[valid] * Bi[valid] / Tot[valid])
  
  # The rate of joining a group for someone not yet part of a group is c/m.
  # The rate of leaving a group is
  join_term <- (c_param / m) * (n - Tot)/n
  
  # The rate of leave a group is B+ or B- depending on the thershold
  frac_red <- ifelse(Tot > 0, Ri / Tot, 0)
  frac_blue <- ifelse(Tot > 0, Bi / Tot, 0)
  
  leaveR_rate <- ifelse(frac_red < T_threshold, beta_plus * Ri, beta_minus * Ri)
  leaveB_rate <- ifelse(frac_blue < T_threshold, beta_plus * Bi, beta_minus * Bi)

  # The rate at which the process leaves the state
  lambda_i <- voter_term + join_term + leaveR_rate + leaveB_rate
  lambda_tot <- sum(lambda_i)
  if (lambda_tot <= 0) break

  # Random group selection
  group_i <- if(m > 1) sample(1:m, 1L, prob = lambda_i / lambda_tot) else 1

  Ri_g <- Ri[group_i]; Bi_g <- Bi[group_i]
  voter_g <- voter_term[group_i]; join_g <- join_term[group_i]
  leaveR_g <- leaveR_rate[group_i]; leaveB_g <- leaveB_rate[group_i]

  rates_vec <- c(ifelse(voter_g>0, voter_g/2, 0),
                 ifelse(voter_g>0, voter_g/2, 0),
                 join_g, leaveR_g, leaveB_g)
  if (sum(rates_vec)<=0) next

  # Apply a random move in Gillespie time

  # Gillespie time
  dt <- rexp(1, lambda_tot)
  t <- t + dt
  if (t >= t_max) break

  # Move step
  # 1: Red to Blue. 2: Blue to Red. 3: Add external to group.
  # 4: Remove internal Red 5: Remove internal Blue
  move <- sample(1:5L, 1L, prob = rates_vec / sum(rates_vec))

  if (move==1 && Ri_g>0) { # Red to Blue
    chosen <- sample(red_members[[group_i]],1)
    opinions[chosen] <- +1
    # update every group where individual belongs
    for (g in groups_of_individual[[chosen]]) {
      idx <- which(red_members[[g]]==chosen)
      if (length(idx)>0) {
        red_members[[g]][idx] <- red_members[[g]][length(red_members[[g]])] # Replaces the chosen individual with the last element
        red_members[[g]] <- red_members[[g]][-length(red_members[[g]])] # delete last
        blue_members[[g]] <- c(blue_members[[g]], chosen)
        Ri[g] <- Ri[g]-1
        Bi[g] <- Bi[g]+1
      }
    }

  } else if (move==2 && Bi_g>0) { # Blue to Red
    chosen <- sample(blue_members[[group_i]],1)
    opinions[chosen] <- -1
    # update every group where individual belongs
    for (g in groups_of_individual[[chosen]]) {
      idx <- which(blue_members[[g]]==chosen)
      if (length(idx)>0) {
        blue_members[[g]][idx] <- blue_members[[g]][length(blue_members[[g]])]
        blue_members[[g]] <- blue_members[[g]][-length(blue_members[[g]])]
        red_members[[g]] <- c(red_members[[g]], chosen)
        Bi[g] <- Bi[g]-1
        Ri[g] <- Ri[g]+1
      }
    }

  } else if (move==3 && length(outsiders[[group_i]])>0) { # Join
    chosen_idx <- sample.int(length(outsiders[[group_i]]),1)
    chosen <- outsiders[[group_i]][chosen_idx]
    outsiders[[group_i]][chosen_idx] <- outsiders[[group_i]][length(outsiders[[group_i]])]
    outsiders[[group_i]] <- outsiders[[group_i]][-length(outsiders[[group_i]])]

    members[[group_i]] <- c(members[[group_i]], chosen)
    groups_of_individual[[chosen]] <- c(groups_of_individual[[chosen]], group_i)
    if (opinions[chosen]==-1) {
      red_members[[group_i]] <- c(red_members[[group_i]], chosen)
      Ri[group_i] <- Ri[group_i]+1
    } else {
      blue_members[[group_i]] <- c(blue_members[[group_i]], chosen)
      Bi[group_i] <- Bi[group_i]+1
    }
    rig_dirty <- TRUE

  } else if (move==4 && Ri_g>0) {
    chosen_idx <- sample.int(Ri_g,1)
    chosen <- red_members[[group_i]][chosen_idx]
    red_members[[group_i]][chosen_idx] <- red_members[[group_i]][length(red_members[[group_i]])]
    red_members[[group_i]] <- red_members[[group_i]][-length(red_members[[group_i]])]

    idx_m <- which(members[[group_i]]==chosen)
    members[[group_i]][idx_m] <- members[[group_i]][length(members[[group_i]])]
    members[[group_i]] <- members[[group_i]][-length(members[[group_i]])]
    groups_of_individual[[chosen]] <- setdiff(groups_of_individual[[chosen]], group_i)
    outsiders[[group_i]] <- c(outsiders[[group_i]], chosen)
    Ri[group_i] <- Ri[group_i]-1
    rig_dirty <- TRUE

  } else if (move==5 && Bi_g>0) {
    chosen_idx <- sample.int(Bi_g,1)
    chosen <- blue_members[[group_i]][chosen_idx]
    blue_members[[group_i]][chosen_idx] <- blue_members[[group_i]][length(blue_members[[group_i]])]
    blue_members[[group_i]] <- blue_members[[group_i]][-length(blue_members[[group_i]])]

    idx_m <- which(members[[group_i]]==chosen)
    members[[group_i]][idx_m] <- members[[group_i]][length(members[[group_i]])]
    members[[group_i]] <- members[[group_i]][-length(members[[group_i]])]
    groups_of_individual[[chosen]] <- setdiff(groups_of_individual[[chosen]], group_i)
    outsiders[[group_i]] <- c(outsiders[[group_i]], chosen)
    Bi[group_i] <- Bi[group_i]-1
    rig_dirty <- TRUE
  }


  # ---- record step ----
  event_counter <- event_counter + 1
  
  if (event_counter %% record_every == 0) {
    opinion_history <- cbind(opinion_history, opinions)
  }
  
  if (event_counter %% record_every == 0) {
    frac_temp <- sapply(levels, function(op) {
      sum(opinions == op)/n
    })
    frac_mat <- cbind(frac_mat, frac_temp)
  }


}

# Final opinion history
opinion_history <- cbind(opinion_history, opinions)

# reconstruct_bipartite <- function(members, n, m) {
#   i <- integer(0); j <- integer(0)
#   for (g in seq_len(m)) {
#     ids <- members[[g]]
#     if (length(ids)>0) { i <- c(i,ids); j <- c(j, rep.int(g,length(ids))) }
#   }
#   sparseMatrix(i=i, j=j, x=1L, dims=c(n,m))
# }
# 
# if (rig_dirty) {
#   B <- reconstruct_bipartite(members, n, m)
#   RIG <- bipartite_to_rig(B)
# }

# list(
#   B0=B0, B=bipartite, RIG=RIG, RIG0=RIG0,
#   #opinions=opinions, #opinions0=opinions0,
#   opinion_history=opinion_history,
#   num_opinions=num_opinions,members=members
# )


#==========================
# --- Output ---
#==========================
#visual_step_multi(RIG0, opinion_history[,1], num_opinions)
#visual_step_multi(RIG, opinion_history[,ncol( opinion_history)], num_opinions)
#visual_bipartite(B0, opinion_history[,1], num_opinions)
#visual_bipartite(B, opinion_history[,ncol( opinion_history)], num_opinions)
visual_step_time(opinion_history, num_opinions)
visual_histo(opinion_history, num_opinions)
par(mfrow =c(1,1))
#heatmapPlot(opinion_history, num_opinions)
visual_histo_pergroup(opinion_history[,1], num_opinions, members0)
visual_histo_pergroup(opinion_history[,ncol( opinion_history)], num_opinions, members)

