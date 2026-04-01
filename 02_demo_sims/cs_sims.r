# Author - Athma
# Last updates - 4/1/2026
# Brief description - Simulates a community of multiple plant species and microbial species

# Inputs - np: number of plant species, nm: number of microbial species, corr_m: correlation matrix, viz: 0/1 to plot

# Outputs - timeseries plot if viz ==1, matrix of covariance, coupling and variance ratio for both plants and microbes, plants only and microbes only

# Notes
## coupling is measured as mean values of the correlation (sub)matrix
## (1/20) scaling standardizes Brownian motion is a slower timescale (400 times = 1/sqrt(400)
## Code translated from MATLAB using ChatGPT 5.2 and checked line-by-line by Athma
## issues - fails for more than 1 plant and 1 microbe

cs_sims <- function(np, nm, corr_m, viz = 0) {
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Install 'mvtnorm': install.packages('mvtnorm')")
  }
  if (!requireNamespace("matrixStats", quietly = TRUE)) {
    stop("Install 'matrixStats': install.packages('matrixStats')")
  }

	# comment if you want replicates
	set.seed(1)

  n <- np + nm
  stopifnot(is.matrix(corr_m), all(dim(corr_m) == c(n, n)))

  # (Optional) make exactly symmetric if corr_m has tiny numerical asymmetry
  corr_m <- (corr_m + t(corr_m)) / 2

  Y <- 10		# total time
  tp <- 0:Y	# timescale for plant dynamics
  tm <- seq(0, Y, by = 1/400)	# timescale for microbial dynamics (roughly daily)

  steps <- length(tm) - 1L  # number of increments
  step_per_unit <- as.integer(round(1 / (tm[2] - tm[1])))  # 400 samples between two plant samples

  # Draw all samples at once (rows correspond to i = 2..length(tm))
  samps <- mvtnorm::rmvnorm(
    n     = steps,
    mean  = rep(0, n),
    sigma = corr_m,
    method = "eigen"   # robust if Sigma is only semi-definite
  )

  xp <- matrix(1, nrow = np, ncol = length(tp))
  xm <- matrix(1, nrow = nm, ncol = length(tm))
  ym <- matrix(0, nrow = nm, ncol = length(tp))

  # Microbes change every tm step
  if (nm > 0) {
    inc_m <- (1/20) * samps[, (np + 1):n, drop = FALSE]          # steps x nm
    xm[, -1] <- 1 + t(matrixStats::colCumsums(inc_m))            # nm x steps
  }

  # Plants update once per integer time (t = 1..Y), using the same samp at that step
  idx <- step_per_unit * (1:Y)                                   # 400, 800, ..., 4000
  if (np > 0) {
    inc_p <- samps[idx, 1:np, drop = FALSE]                       # Y x np
    xp[, -1] <- 1 + t(matrixStats::colCumsums(inc_p))             # np x Y
  }

  # Record microbes at those integer times (after xm update)
  if (nm > 0) {
    ym[, -1] <- xm[, idx + 1L, drop = FALSE]                      # +1 because xm includes time 0
  }

  if (viz == 1) {
		dev.new()
		plot(tp,xp, type = "o", pch = 1, lty = 1,xlab = "Time", ylab = "Variables",ylim=c(min(xp,xm),max(xp,xm)))
		lines(tp,xp)
		lines(tm,xm)
  }

	# cov - covariance, vr - variance ration, cpl - coupling
	# pm - plants and microbes, p - plants only, m - microbes only

  cov_pm <- stats::cov(t(rbind(xp, ym)))
  vr_pm  <- (sum(cov_pm) - sum(diag(cov_pm))) / sum(diag(cov_pm))
  cpl_pm <- mean(abs(corr_m))

  cov_p <- stats::cov(t(xp))
  vr_p  <- (sum(cov_p) - sum(diag(cov_p))) / sum(diag(cov_p))
  cpl_p <- mean(abs(corr_m[1:np, 1:np, drop = FALSE]))

  cov_m <- stats::cov(t(xm))
  vr_m  <- (sum(cov_m) - sum(diag(cov_m))) / sum(diag(cov_m))
  cpl_m <- mean(abs(corr_m[(np + 1):n, (np + 1):n, drop = FALSE]))

	print(cov_pm)

  out <- rbind(
 #   cov = c(cov_pm, cov_p, cov_m),
		cpl = c(cpl_pm, cpl_p, cpl_m),
    vr  = c(vr_pm,  vr_p,  vr_m)
  )
  colnames(out) <- c("pm", "p", "m")
  out
}
