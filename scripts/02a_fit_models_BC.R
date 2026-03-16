################################################################################
# 02a_fit_models_BC.R — BC-Only Dynamic Occupancy Models (all species)
#
# Loops through all species and fits each one.
# Input:  data/processed/BC/
# Output: outputs/BC/fits/
#
# Covariates:
#   ψ (initial occupancy): Elevation + Water Nearby + Region (random effect)
#   γ (colonization):      Distance to Road + Distance to Harvest
#   φ (persistence):       Distance to Road + Distance to Harvest
#   p (detection):         Clutter + Temperature + Julian Date
################################################################################

library(jagsUI)

# ============================================================================
# SETTINGS
# ============================================================================

SPECIES_LIST <- c("ANPA", "COTO", "EPFU", "EUMA", "LABO", "LACI", "LANO",
                  "MYCA", "MYCI", "MYEV", "MYLU", "MYSE", "MYTH", "MYVO",
                  "MYYU", "PAHE", "TABR")

INPUT_DIR  <- "data/processed/BC"
OUTPUT_DIR <- "outputs/BC/fits"

N_ITER <- 30000; N_BURNIN <- 15000; N_THIN <- 5; N_CHAINS <- 3

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# JAGS MODEL CODE (same for all species)
# ============================================================================

model_code <- "
model {
  beta.psi[1] ~ dnorm(0, 0.368)
  beta.psi[2] ~ dnorm(0, 0.368)
  beta.psi[3] ~ dnorm(0, 0.368)

  sigma.region ~ dunif(0, 5)
  tau.region <- 1 / (sigma.region * sigma.region)
  for (r in 1:nregion) { u.region[r] ~ dnorm(0, tau.region) }

  for (t in 1:(nyear-1)) { alpha.gamma[t] ~ dnorm(0, 0.368) }
  beta.gamma[1] ~ dnorm(0, 0.368)
  beta.gamma[2] ~ dnorm(0, 0.368)

  for (t in 1:(nyear-1)) { alpha.phi[t] ~ dnorm(0, 0.368) }
  beta.phi[1] ~ dnorm(0, 0.368)
  beta.phi[2] ~ dnorm(0, 0.368)

  for (k in 1:4) { beta.p[k] ~ dnorm(0, 0.368) }

  for (i in 1:nsite) {
    logit(psi[i]) <- beta.psi[1] + beta.psi[2]*elevation[i] +
                     beta.psi[3]*water[i] + u.region[region[i]]
    z[i, 1] ~ dbern(psi[i])
    for (t in 2:nyear) {
      logit(gamma[i, t-1]) <- alpha.gamma[t-1] +
                              beta.gamma[1]*dist_road[i] +
                              beta.gamma[2]*dist_harvest[i]
      logit(phi[i, t-1])   <- alpha.phi[t-1] +
                              beta.phi[1]*dist_road[i] +
                              beta.phi[2]*dist_harvest[i]
      muZ[i, t] <- z[i, t-1]*phi[i, t-1] + (1 - z[i, t-1])*gamma[i, t-1]
      z[i, t] ~ dbern(muZ[i, t])
    }
  }

  for (i in 1:nsite) {
    for (t in 1:nyear) {
      for (j in 1:J[i, t]) {
        logit(p[i, t, j]) <- beta.p[1] +
                             beta.p[2]*clutter[i, t, nsurv[i, t, j]] +
                             beta.p[3]*temp[i, t, nsurv[i, t, j]] +
                             beta.p[4]*julian[i, t, nsurv[i, t, j]]
        y[i, t, nsurv[i, t, j]] ~ dbern(z[i, t] * p[i, t, j])
      }
    }
  }

  for (t in 1:nyear) { psi.fs[t] <- mean(z[, t]); n.occ[t] <- sum(z[, t]) }
  for (t in 1:(nyear-1)) {
    gamma.avg[t] <- ilogit(alpha.gamma[t])
    phi.avg[t]   <- ilogit(alpha.phi[t])
    ext.avg[t]   <- 1 - phi.avg[t]
  }
  mean.gamma <- mean(gamma.avg[])
  mean.phi   <- mean(phi.avg[])
  mean.ext   <- 1 - mean.phi
  lambda     <- psi.fs[nyear] / max(psi.fs[1], 0.001)
  for (t in 1:nyear) { year.centered[t] <- t - (nyear + 1) / 2 }
  trend.slope <- sum(year.centered[] * psi.fs[]) /
                 sum(year.centered[] * year.centered[])
  mean.p <- ilogit(beta.p[1])
}
"

model_file <- file.path(OUTPUT_DIR, "model_dynocc_BC.txt")
writeLines(model_code, model_file)

params <- c("beta.psi", "beta.gamma", "beta.phi", "beta.p",
            "alpha.gamma", "alpha.phi", "gamma.avg", "phi.avg", "ext.avg",
            "sigma.region", "u.region", "psi.fs", "n.occ",
            "mean.gamma", "mean.phi", "mean.ext", "mean.p",
            "lambda", "trend.slope", "z")

# ============================================================================
# LOOP THROUGH SPECIES
# ============================================================================

for (sp in SPECIES_LIST) {

  cat("\n===================================================\n")
  cat("BC Model:", sp, "\n")
  cat("===================================================\n")

  # Load data
  data_file <- file.path(INPUT_DIR, paste0("QUAD_", sp, "_data.rds"))
  if (!file.exists(data_file)) { cat("  Data not found, skipping\n"); next }
  d <- readRDS(data_file)
  nyear <- d$nyear

  cat("Sites:", d$nsite, "| Years:", nyear,
      "(", min(d$years), "-", max(d$years), ")\n")
  cat("Detections:", sum(d$y == 1, na.rm = TRUE),
      "| Surveys:", sum(!is.na(d$y)), "\n")

  # Prepare JAGS input
  temp_clean    <- d$temp;    temp_clean[is.na(temp_clean)]       <- 0
  clutter_clean <- d$clutter; clutter_clean[is.na(clutter_clean)] <- 0
  julian_clean  <- d$julian;  julian_clean[is.na(julian_clean)]   <- 0

  jags_data <- list(
    y = d$y, J = d$J, nsurv = d$nsurv,
    nsite = d$nsite, nyear = nyear, nregion = d$nregion,
    elevation = d$elevation, water = d$water,
    dist_road = d$dist_road, dist_harvest = d$dist_harvest,
    region = d$region,
    temp = temp_clean, clutter = clutter_clean, julian = julian_clean
  )

  init_z <- apply(d$y, 1:2, function(x) {
    if (all(is.na(x))) return(1)
    ifelse(any(x == 1, na.rm = TRUE), 1, rbinom(1, 1, 0.5))
  })

  inits <- function() list(
    z = init_z, beta.psi = rnorm(3, 0, 0.5),
    beta.gamma = rnorm(2, 0, 0.5), beta.phi = rnorm(2, 0, 0.5),
    beta.p = rnorm(4, 0, 0.5),
    alpha.gamma = rnorm(nyear - 1, -1, 0.5),
    alpha.phi = rnorm(nyear - 1, 1, 0.5),
    sigma.region = runif(1, 0.1, 1)
  )

  # Fit
  cat("Fitting (", N_ITER, " iter, ", N_CHAINS, " chains)...\n", sep = "")
  start_time <- Sys.time()

  fit <- jags(data = jags_data, inits = inits, parameters.to.save = params,
              model.file = model_file,
              n.chains = N_CHAINS, n.iter = N_ITER, n.burnin = N_BURNIN,
              n.thin = N_THIN, parallel = TRUE)

  elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 1)
  cat("Done in", elapsed, "min\n")

  # Save
  out_tag <- paste0("QUAD_", sp)
  saveRDS(fit, file.path(OUTPUT_DIR, paste0(out_tag, "_fit.rds")))

  summ <- fit$summary
  summ_df <- as.data.frame(summ); summ_df$parameter <- rownames(summ)
  write.csv(summ_df, file.path(OUTPUT_DIR, paste0(out_tag, "_summary.csv")),
            row.names = FALSE)

  # Print summary
  cat("\nRESULTS:", sp, "(BC)\n")
  cat("Convergence: max Rhat =", round(max(summ[,"Rhat"], na.rm=TRUE), 3),
      "| min n.eff =", round(min(summ[,"n.eff"], na.rm=TRUE)), "\n")

  lam <- summ["lambda", "mean"]
  lam_lo <- summ["lambda", "2.5%"]; lam_hi <- summ["lambda", "97.5%"]
  cat(sprintf("  lambda = %.2f (%.2f, %.2f)", lam, lam_lo, lam_hi))
  if (lam_hi < 1) cat(" -> DECLINE\n") else
    if (lam_lo > 1) cat(" -> INCREASE\n") else
      cat(" -> No clear trend\n")
}

cat("\n\nAll BC models complete. Output:", OUTPUT_DIR, "\n")
