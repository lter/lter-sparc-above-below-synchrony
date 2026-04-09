library(mvtnorm)

cs_sims <- function(np, nm, intm, viz) {
  n <- np + nm
  Y <- 10
  
  tp <- seq(0, Y, by = 1)
  tm <- seq(0, Y, by = 1/400)
  
  xp <- matrix(1, nrow = np, ncol = length(tp))
  xm <- matrix(1, nrow = nm, ncol = length(tm))
  ym <- matrix(0, nrow = nm, ncol = length(tp))
  
  tpi <- 1
  
  # Simulation Loop
  for (i in 2:length(tm)) {
    # Generate random sample
    samp <- rmvnorm(1, mean = rep(0, n), sigma = intm)
    
    # Update microbes
    xm[, i] <- xm[, i - 1] + (1/20) * samp[(np + 1):n]
    
    # Check if current tm matches the next tp index
    # Using a small epsilon to avoid floating point issues
    if (tpi < length(tp) && abs(tm[i] - tp[tpi + 1]) < 1e-10) {
      tpi <- tpi + 1
      xp[, tpi] <- xp[, tpi - 1] + samp[1:np]
      ym[, tpi] <- xm[, i]
    }
  }
  
  # Visualization
	if (viz == 1) {
		dev.new()
    plot(tp, xp[1, ], type = "b", col = "blue", ylim = range(c(xp, xm)), 
         xlab = "Time", ylab = "Variables", main = "Simulation Viz")
    for(r in 2:np) lines(tp, xp[r, ], type = "b", col = "blue")
    for(r in 1:nm) lines(tm, xm[r, ], col = "gray")
  }
  
  # Statistics Calculations
  calc_vr <- function(cv) {
    (sum(cv) - sum(diag(cv))) / sum(diag(cv))
  }
  
  # Combined (Plants + Microbes)
  cov_pm <- cov(t(rbind(xp, ym)))
  vr_pm  <- calc_vr(cov_pm)
  cpl_pm <- mean(abs(intm))
  
  # Plants only
  cov_p  <- cov(t(xp))
  vr_p   <- calc_vr(cov_p)
  cpl_p  <- mean(abs(intm[1:np, 1:np]))
  
  # Microbes only
  cov_m  <- cov(t(xm))
  vr_m   <- calc_vr(cov_m)
  cpl_m  <- mean(abs(intm[(np + 1):n, (np + 1):n]))
  
  # Return matrix matches MATLAB output structure
  return(matrix(c(cpl_pm, vr_pm, cpl_p, vr_p, cpl_m, vr_m), nrow = 2, ncol = 3))
}
