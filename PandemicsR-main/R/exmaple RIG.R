
#Graph 1: Two clusters, opposing opinions connected
library(Matrix)

n <- 15; m <- 3
B <- Matrix(0, n, m, sparse = TRUE)

# Cluster 1: nodes 1-5 in group 1
B[1:5, 1] <- 1
# Cluster 2: nodes 6-10 in group 2
B[6:10, 2] <- 1
# Cluster 3: nodes 11-15 in group 3
B[11:15, 3] <- 1

# Cross connections to allow interaction
B[5,2] <- 1
B[10,3] <- 1
B[1,3] <- 1

# Generate RIG
RIG <- (B %*% t(B)) > 0
diag(RIG) <- 0

# Assign opinions
opinions <- c(rep(1,5), rep(-1,5), rep(1,5))


#Alternating opinions with cross-group edges
B <- Matrix(0, n, m, sparse = TRUE)

# Spread nodes across groups
B[1:5,1] <- 1
B[6:10,2] <- 1
B[11:15,3] <- 1

# Add cross-group edges for interaction
B[c(3,8,13), c(2,3,1)] <- 1

# Generate RIG
RIG <- (B %*% t(B)) > 0
diag(RIG) <- 0

# Assign alternating opinions
opinions <- rep(c(1,-1,1), each = 5)
