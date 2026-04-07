library(mvtnorm)

cs_sims <- function(np, nm, intm, viz) {
  n <- np + nm # total number of species
  Y <- 10 # number of years
  
  tp <- seq(0, Y, by = 1) # time points for plants
  tm <- seq(0, Y, by = 1/400) # time points for microbes
  
  # initialize abundance matrices
  xp <- matrix(1, nrow = np, ncol = length(tp)) # plants (updated yearly)
  xm <- matrix(1, nrow = nm, ncol = length(tm)) # microbes (updated subyearly)
  ym <- matrix(0, nrow = nm, ncol = length(tp)) # microbes (tracked yearly)
  
  tpi <- 1 # initialize plant time step tracker
  
  # Simulation Loop
  for (i in 2:length(tm)) {
    # Generate random abundance fluctuations
    samp <- rmvnorm(1, mean = rep(0, n), sigma = intm)
    
    # Update microbes
    xm[, i] <- xm[, i - 1] + (1/20) * samp[(np + 1):n]
    
    # Check if current tm matches the next tp index
    # Using a small epsilon to avoid floating point issues
    if (tpi < length(tp) && abs(tm[i] - tp[tpi + 1]) < 1e-10) {
      tpi <- tpi + 1 # update plant time step tracker
      xp[, tpi] <- xp[, tpi - 1] + samp[1:np] # update plant abundance
      ym[, tpi] <- xm[, i] # record annual microbe abundance
    }
  }
  
  # Timeseries plots for plants and microbes
  if (viz == 1) {
    plot(tp, xp[1, ], type = "b", col = "blue", ylim = range(c(xp, xm)), 
         xlab = "Time", ylab = "Variables", main = "Simulation Viz")
    for(r in 2:np) lines(tp, xp[r, ], type = "b", col = "blue")
    for(r in 1:nm) lines(tm, xm[r, ], col = "gray", alpha = 0.5)
  }
  
  # Statistics Calculations
  calc_vr <- function(cv) { # Variance Ratio
    (sum(cv) - sum(diag(cv))) / sum(diag(cv))
  }
  
  # Combined (Plants + Microbes)
  cov_pm <- cov(t(rbind(xp, ym))) # Timeseries Covariance
  vr_pm  <- calc_vr(cov_pm) # Variance Ratio
  cpl_pm <- mean(abs(intm)) # Degree of Coupling
  
  # Plants only
  cov_p  <- cov(t(xp)) # Timeseries Covariance
  vr_p   <- calc_vr(cov_p) # Variance Ratio
  cpl_p  <- mean(abs(intm[1:np, 1:np])) # Degree of Coupling
  
  # Microbes only
  cov_m  <- cov(t(xm)) # Timeseries Covariance
  vr_m   <- calc_vr(cov_m) # Variance Ratio
  cpl_m  <- mean(abs(intm[(np + 1):n, (np + 1):n])) # Degree of Coupling
  
  # Return matrix matches MATLAB output structure
  return(matrix(c(cpl_pm, vr_pm, cpl_p, vr_p, cpl_m, vr_m), nrow = 2, ncol = 3))
}
