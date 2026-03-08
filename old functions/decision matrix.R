# Define the decision matrix as a data.frame
decision_df <- data.frame(
  DR = c(0.0, 0.1, 0.05, 0.05),
  LR = c(0.2, 0.0, 0.15, 0.1),
  LB = c(0.15, 0.15, 0.0, 0.15),
  DB = c(0.05, 0.05, 0.1, 0.0),
  row.names = c("DR", "LR", "LB", "DB"),
  check.names = FALSE
)
rownames(decision_df)
# Inspect
decision_df <- decision_df*0

decision_df[2,1] <- 0.1 # light to dark same color c1
decision_df[3,4] <- 0.1 # light to dark same color c1

decision_df[3,1] <- 0.2 # light to dark diff color c2
decision_df[2,4] <- 0.2 # light to dark diff color c2

decision_df[2,3] <- 0.3 # light to light diff color c3
decision_df[3,2] <- 0.3 # light to light diff color c3

decision_df[3,4] <- 0.4 # light to light diff color c4
decision_df[2,1] <- 0.4 # light to light diff color c4

decision_df
