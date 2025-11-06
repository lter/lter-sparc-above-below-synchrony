#GPT generated from my matlab code. The comments are good but the code is bad!

# Clear workspace
rm(list = ls())
graphics.off()

set.seed(1)

dt <- 0.1      # resolution
T <- 100       # length
N <- 2        # number of variables

t <- seq(0, T, by = dt)
nt <- length(t)

ns1 <- 0.5
ns2 <- 1.5

# Generate random matrices
x <- sin(0.1 * t + pi * ns2 * matrix(runif(N), nrow = N, ncol = 1)) + ns1 * matrix(runif(N * nt), nrow = N, ncol = nt)

# Covariance and correlation (note: R's cov uses columns as variables)
cm <- cov(t(x))
corr_mat <- cor(t(x))

# Variance ratio calculation
vr <- (sum(cm) - sum(diag(cm))) / sum(diag(cm))
print(vr)

# Downsample: take every 10th column
idx <- seq(1, nt, by = 10)
y <- x[, idx]

# Covariance and variance ratio for downsampled data
cm <- cov(t(y))
vr <- (sum(cm) - sum(diag(cm))) / sum(diag(cm))
print(vr)

# Plotting
# The *- marker is achieved by type="b" and pch=8 or similar
matplot(t, t(x), type = "b", pch = 8, lty = 1, xlab = "t", ylab = "x", main = "All Variables")
matplot(t[idx], t(y), type = "b", pch = 8, lty = 1, xlab = "t (downsampled)", ylab = "y", main = "Downsampled Variables")
