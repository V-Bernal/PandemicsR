#==========================
# Libraries
#==========================

library(igraph)
library(Matrix)
library(PandemicsR)
library(zip)


#==========================
# Load project functions
#==========================

# If your functions are in separate files:
# source("simulation.R")
# source("visualization.R")


#==========================
# Parameters
#==========================

params <- list(

  # Network
  n = 15,
  m = 3,
  t_max = 100,
  lambda = 1,


  # Schelling
  runSchelling = FALSE,
  c_param = 0.4,
  beta_plus = 0.5,
  beta_minus = 0.2,
  T_threshold = 0.3,


  # Voter
  runVoter = TRUE,
  gamma = 5,
  num_opinions = 4,

  alpha = 0.1,
  alpha_deradicalization = 0,
  alpha0 = 0,


  # Epidemic
  runEpidemic = FALSE,
  I0 = 0.01,
  gamma_epi = 0.3,
  beta_red = 0.5,
  beta_blue = 0.2

)



#==========================
# Run simulation
#==========================

start_time <- Sys.time()

sim <- run_simulation(params)

end_time <- Sys.time()


cat(
  "Simulation stopped:",
  sim$stop_reason,
  "\n"
)

cat(
  "Final Gillespie time:",
  sim$final_time,
  "\n"
)

cat(
  "Number of events:",
  sim$event_counter,
  "\n"
)

cat(
  "Computation time:",
  difftime(end_time,start_time,units="secs"),
  "seconds\n"
)



#==========================
# Generate plots
#==========================

dir.create(
  "plots",
  showWarnings = FALSE
)


# Opinion dynamics

png(
  "plots/opinion_fraction.png",
  width=1400,
  height=1000,
  res=150
)

visual_step_time(
  sim$frac_mat,
  sim$num_opinions
)

dev.off()



# Histogram

png(
  "plots/opinion_histogram.png",
  width=1400,
  height=1000,
  res=150
)

visual_histo(
  sim$opinion_history,
  sim$num_opinions
)

dev.off()



# Group histogram

png(
  "plots/group_histogram_initial.png",
  width=1400,
  height=1000,
  res=150
)

visual_histo_pergroup(
  sim$opinion_history[,1],
  sim$num_opinions,
  sim$members0
)

dev.off()



png(
  "plots/group_histogram_final.png",
  width=1400,
  height=1000,
  res=150
)

visual_histo_pergroup(
  sim$opinion_history[,ncol(sim$opinion_history)],
  sim$num_opinions,
  sim$members
)

dev.off()



#==========================
# Graph visualization
#==========================


png(
  "plots/final_RIG.png",
  width=1400,
  height=1000,
  res=150
)

visual_step_multi(
  sim$RIG,
  sim$opinion_history[,ncol(sim$opinion_history)],
  sim$num_opinions
)

dev.off()



png(
  "plots/final_bipartite.png",
  width=1400,
  height=1000,
  res=150
)

visual_bipartite(
  sim$B,
  sim$opinion_history[,ncol(sim$opinion_history)],
  sim$num_opinions
)

dev.off()



#==========================
# Epidemic plots
#==========================

if(sim$runEpidemic){

  png(
    "plots/SIR.png",
    width=1400,
    height=1000,
    res=150
  )

  visual_step_time_SIR(
    sim$SIR_df
  )

  dev.off()

}



#==========================
# Save simulation object
#==========================

saveRDS(
  sim,
  "simulation_results.rds"
)
