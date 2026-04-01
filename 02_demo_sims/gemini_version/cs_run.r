source("cs_sims.r")

# Setup
set.seed(1) # Equivalent to rng(1)

np <- 5
nm <- 10
n  <- np + nm

start_time <- Sys.time()

intm <- matrix(rnorm(n * n, mean = 0, sd =  0.01), nrow = n)
diag(intm) <- 1
intm <- intm %*% t(intm) # Ensure symmetry/Positive Definiteness

out <- cs_sims(np, nm, intm, 1)

print(out)

print(Sys.time() - start_time)
