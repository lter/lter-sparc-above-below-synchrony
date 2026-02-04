cs_sims <- function(np, nm, intm, viz = 0) {
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Install 'mvtnorm': install.packages('mvtnorm')")
  }
  if (!requireNamespace("matrixStats", quietly = TRUE)) {
    stop("Install 'matrixStats': install.packages('matrixStats')")
  }

	# comment if you want replicates
	set.seed(1)

  n <- np + nm
  stopifnot(is.matrix(intm), all(dim(intm) == c(n, n)))

  # (Optional) make exactly symmetric if intm has tiny numerical asymmetry
  intm <- (intm + t(intm)) / 2

  Y <- 10
  tp <- 0:Y
  tm <- seq(0, Y, by = 1/400)

  steps <- length(tm) - 1L                 # number of increments
  step_per_unit <- as.integer(round(1 / (tm[2] - tm[1])))  # 400

  # Draw all samples at once (rows correspond to i = 2..length(tm))
  samps <- mvtnorm::rmvnorm(
    n     = steps,
    mean  = rep(0, n),
    sigma = intm,
    method = "eigen"   # robust if Sigma is only semi-definite
  )

  xp <- matrix(1, nrow = np, ncol = length(tp))
  xm <- matrix(1, nrow = nm, ncol = length(tm))
  ym <- matrix(0, nrow = nm, ncol = length(tp))

  # Microbes evolve every tm step
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
    graphics::matplot(tp, t(xp), type = "o", pch = 1, lty = 1,
                      xlab = "Time", ylab = "Variables")
    graphics::matlines(tm, t(xm), lty = 1)
  }

  cov_pm <- stats::cov(t(rbind(xp, ym)))
  vr_pm  <- (sum(cov_pm) - sum(diag(cov_pm))) / sum(diag(cov_pm))
  cpl_pm <- mean(abs(intm))

  cov_p <- stats::cov(t(xp))
  vr_p  <- (sum(cov_p) - sum(diag(cov_p))) / sum(diag(cov_p))
  cpl_p <- mean(abs(intm[1:np, 1:np, drop = FALSE]))

  cov_m <- stats::cov(t(xm))
  vr_m  <- (sum(cov_m) - sum(diag(cov_m))) / sum(diag(cov_m))
  cpl_m <- mean(abs(intm[(np + 1):n, (np + 1):n, drop = FALSE]))

  out <- rbind(
    cpl = c(cpl_pm, cpl_p, cpl_m),
    vr  = c(vr_pm,  vr_p,  vr_m)
  )
  colnames(out) <- c("pm", "p", "m")
  out
}
