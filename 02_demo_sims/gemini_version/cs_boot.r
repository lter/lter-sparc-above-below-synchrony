source("cs_sims.r")

# Setup
set.seed(1) # Equivalent to rng(1)

np <- 5 # number of plants
nm <- 10 # number of microbes
n  <- np + nm # total number of species
nc <- 20 # number of communities
bootn <- 10 # number of replicates

# Initialize 4D array: [row, col, community, replicate]
summ <- array(0, dim = c(2, 3, nc, bootn))

start_time <- Sys.time() # to track performance

for (i in 1:nc) {
  # Create interaction matrix
  intm <- matrix(rnorm(n * n, mean = 0, sd = i * 0.01), nrow = n)
  diag(intm) <- 1
  intm <- intm %*% t(intm) # Ensure symmetry/Positive Definiteness
  
  for (k in 1:bootn) {
    # Visualize only on the first community, last bootstrap replicate
    do_viz <- if (i == 1 && k == bootn) 1 else 0
    summ[, , i, k] <- cs_sims(np, nm, intm, do_viz)
  }
  print(paste("Community:", i))
}

print(Sys.time() - start_time) # how long did the simulations take? 

# --- Plotting Results ---

# Flatten the dimensions for plotting
# MATLAB: reshape(squeeze(summ(1,1,:,:)), [nc*bootn 1])
coupling <- as.vector(summ[1, 1, , ])
vr_pm    <- as.vector(summ[2, 1, , ])
vr_p     <- as.vector(summ[2, 2, , ])
vr_m     <- as.vector(summ[2, 3, , ])

plot(coupling, vr_pm, pch = 8, col = "black", xlab = "Coupling", ylab = "Variance ratio", ylim = range(c(vr_pm, vr_p, vr_m)))
points(coupling, vr_p, pch = 8, col = "red")
points(coupling, vr_m, pch = 8, col = "blue")
legend("topright", legend = c("pm", "p", "m"), col = c("black", "red", "blue"), pch = 8)
