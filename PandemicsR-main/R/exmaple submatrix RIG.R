# Example recompute RIG sub matrix

library(Matrix)

# Small bipartite matrix (5 individuals × 3 groups)
bipartite <- Matrix(
  c(1,0,1,  # node 1 in groups 1,3
    0,1,1,  # node 2 in groups 2,3
    1,0,0,  # node 3 in group 1
    0,1,0,  # node 4 in group 2
    0,0,1), # node 5 in group 3
  nrow = 5, ncol = 3, byrow = TRUE, sparse = TRUE)

colnames(bipartite) <- paste0("G",1:3)
rownames(bipartite) <- paste0("I",1:5)

# Compute initial RIG (projection)
RIG_full <- (bipartite %*% t(bipartite)) > 0
diag(RIG_full) <- 0
RIG_full <- Matrix(RIG_full, sparse = TRUE)
RIG_full

# Copy original for local update
bip_new <- bipartite
bip_new[2,3]
bip_new[2,3] <- 0  # node 2 leaves group 3

bipartite
bip_new

i <- 2
neighbors_of_i <- which(RIG_full[i, ] == 1)
affected <- sort(unique(c(i, neighbors_of_i)))
affected
# expected: nodes 1,2,4,5 (since 2 was linked to 1,4,5)

# compare computations
sub_RIG <- (bip_new[affected, ] %*% t(bip_new[affected, ])) > 0
diag(sub_RIG) <- 0
RIG_updated <- RIG_full
RIG_updated
RIG_updated[affected, affected] <- Matrix(sub_RIG, sparse = TRUE)

RIG_full_recomputed <- (bip_new %*% t(bip_new)) > 0
diag(RIG_full_recomputed) <- 0
RIG_full_recomputed <- Matrix(RIG_full_recomputed, sparse = TRUE)

all.equal(RIG_updated, RIG_full_recomputed)

as.matrix(RIG_full_recomputed)
as.matrix(RIG_updated)

