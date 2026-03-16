################################################################################
# BC-Only Dynamic Occupancy Model
#
# Fits per-species JAGS model with dist_road + dist_harvest on gamma/phi.
#
# USAGE: Set SPECIES, then source() or run interactively.
#
# Input:  data/processed/BC/
# Output: outputs/BC/fits/
################################################################################

library(jagsUI)

# ============================================================================
# SETTINGS
# ============================================================================

#SPECIES <- "EPFU"   # change per species

SPECIES <- c("ANPA", "COTO", "EPFU", "EUMA", "LABO", "LACI", "LANO",
         "MYCA", "MYCI", "MYEV", "MYLU", "MYSE", "MYTH", "MYVO",
         "MYYU", "PAHE", "TABR")

INPUT_DIR  <- "data/processed/BC"
OUTPUT_DIR <- "outputs/BC/fits"

N_ITER <- 30000; N_BURNIN <- 15000; N_THIN <- 5; N_CHAINS <- 3

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# LOAD DATA
# ============================================================================

# loop through each species
for (i in SPECIES) {
d <- readRDS(file.path(INPUT_DIR, paste0("QUAD_", SPECIES, "_data.rds")))
nyear <- d$nyear

cat("=== BC Model:", SPECIES, "===\n")
cat("Sites:", d$nsite, "| Years:", nyear,
    "(", min(d$years), "-", max(d$years), ")\n")
cat("Regions:", d$nregion, "-", paste(d$regions, collapse = ", "), "\n")
cat("Detections:", sum(d$y == 1, na.rm = TRUE),
    "| Surveys:", sum(!is.na(d$y)), "\n\n")

# ============================================================================
# JAGS MODEL — BC (with dist_road + dist_harvest)
# ============================================================================

model_code <- "
model {
  # Priors
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

  # State process
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

  # Observation process
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

  # Derived
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

# ============================================================================
# PREPARE JAGS INPUT
# ============================================================================

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

params <- c("beta.psi", "beta.gamma", "beta.phi", "beta.p",
            "alpha.gamma", "alpha.phi", "gamma.avg", "phi.avg", "ext.avg",
            "sigma.region", "u.region", "psi.fs", "n.occ",
            "mean.gamma", "mean.phi", "mean.ext", "mean.p",
            "lambda", "trend.slope")

# ============================================================================
# FIT
# ============================================================================

cat("Fitting (", N_ITER, " iter, ", N_CHAINS, " chains)...\n", sep = "")
start_time <- Sys.time()

fit <- jags(data = jags_data, inits = inits, parameters.to.save = params,
            model.file = model_file,
            n.chains = N_CHAINS, n.iter = N_ITER, n.burnin = N_BURNIN,
            n.thin = N_THIN, parallel = TRUE)

cat("Done in", round(difftime(Sys.time(), start_time, units = "mins"), 1), "min\n\n")

# ============================================================================
# SAVE & PRINT
# ============================================================================

out_tag <- paste0("QUAD_", SPECIES)
saveRDS(fit, file.path(OUTPUT_DIR, paste0(out_tag, "_fit.rds")))

summ <- fit$summary
summ_df <- as.data.frame(summ); summ_df$parameter <- rownames(summ)
write.csv(summ_df, file.path(OUTPUT_DIR, paste0(out_tag, "_summary.csv")), row.names = FALSE)

sig <- function(lo, hi) ifelse(lo > 0 | hi < 0, " *", "")
pr <- function(label, param) {
  cat(sprintf("  %-20s %7.3f (%7.3f, %7.3f)%s\n",
              label, summ[param, "mean"], summ[param, "2.5%"], summ[param, "97.5%"],
              sig(summ[param, "2.5%"], summ[param, "97.5%"])))
}

cat("RESULTS:", SPECIES, "(BC)\n")
cat("Convergence: max Rhat =", round(max(summ[,"Rhat"], na.rm=TRUE), 3),
    "| min n.eff =", round(min(summ[,"n.eff"], na.rm=TRUE)), "\n\n")

cat("Initial Occupancy (psi):\n")
pr("Intercept", "beta.psi[1]"); pr("Elevation", "beta.psi[2]")
pr("Water", "beta.psi[3]"); pr("Region sigma", "sigma.region")

cat("\nColonization (gamma):\n")
pr("Dist to Road", "beta.gamma[1]"); pr("Dist to Harvest", "beta.gamma[2]")
for (t in 1:(nyear-1))
  pr(paste0(d$years[t], "->", d$years[t+1]), paste0("gamma.avg[", t, "]"))

cat("\nPersistence (phi):\n")
pr("Dist to Road", "beta.phi[1]"); pr("Dist to Harvest", "beta.phi[2]")
for (t in 1:(nyear-1))
  pr(paste0(d$years[t], "->", d$years[t+1]), paste0("phi.avg[", t, "]"))

cat("\nDetection (p):\n")
pr("Intercept", "beta.p[1]"); pr("Clutter", "beta.p[2]")
pr("Temperature", "beta.p[3]"); pr("Julian", "beta.p[4]")

cat("\nDerived:\n")
pr("Mean colonization", "mean.gamma"); pr("Mean persistence", "mean.phi")
pr("Mean extinction", "mean.ext"); pr("Mean detection", "mean.p")
pr("lambda", "lambda"); pr("Trend slope", "trend.slope")

lam <- summ["lambda", ]
if (lam["97.5%"] < 1) cat("\n  -> DECLINE\n") else
  if (lam["2.5%"] > 1) cat("\n  -> INCREASE\n") else
    cat("\n  -> No clear trend\n")

cat("\nSaved:", file.path(OUTPUT_DIR, paste0(out_tag, "_fit.rds")), "\n")



################################################################################
# BC + Alaska Dynamic Occupancy Model
#
# Fits per-species JAGS model with dist_road only on gamma/phi.
# No dist_harvest (unavailable for Alaska).
#
# NOTE: With only 1 distance covariate, JAGS names the parameters
# "beta.gamma" and "beta.phi" (no [1] index). The print section
# below accounts for this.
#
# USAGE: Set SPECIES, then source() or run interactively.
#
# Input:  data/processed/BCAK/
# Output: outputs/BCAK/fits/
################################################################################

library(jagsUI)

# ============================================================================
# SETTINGS
# ============================================================================

SPECIES <- "EPFU"   # change per species

INPUT_DIR  <- "data/processed/BCAK"
OUTPUT_DIR <- "outputs/BCAK/fits"

N_ITER <- 30000; N_BURNIN <- 15000; N_THIN <- 5; N_CHAINS <- 3

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# LOAD DATA
# ============================================================================

d <- readRDS(file.path(INPUT_DIR, paste0("QUAD_", SPECIES, "_data.rds")))
nyear <- d$nyear

cat("=== BCAK Model:", SPECIES, "===\n")
cat("Sites:", d$nsite, "| Years:", nyear,
    "(", min(d$years), "-", max(d$years), ")\n")
cat("Regions:", d$nregion, "-", paste(d$regions, collapse = ", "), "\n")
cat("Detections:", sum(d$y == 1, na.rm = TRUE),
    "| Surveys:", sum(!is.na(d$y)), "\n\n")

# ============================================================================
# JAGS MODEL — BCAK (dist_road only, no harvest)
# ============================================================================

model_code <- "
model {
  # Priors
  beta.psi[1] ~ dnorm(0, 0.368)
  beta.psi[2] ~ dnorm(0, 0.368)
  beta.psi[3] ~ dnorm(0, 0.368)

  sigma.region ~ dunif(0, 5)
  tau.region <- 1 / (sigma.region * sigma.region)
  for (r in 1:nregion) { u.region[r] ~ dnorm(0, tau.region) }

  for (t in 1:(nyear-1)) { alpha.gamma[t] ~ dnorm(0, 0.368) }
  beta.gamma ~ dnorm(0, 0.368)

  for (t in 1:(nyear-1)) { alpha.phi[t] ~ dnorm(0, 0.368) }
  beta.phi ~ dnorm(0, 0.368)

  for (k in 1:4) { beta.p[k] ~ dnorm(0, 0.368) }

  # State process
  for (i in 1:nsite) {
    logit(psi[i]) <- beta.psi[1] + beta.psi[2]*elevation[i] +
                     beta.psi[3]*water[i] + u.region[region[i]]
    z[i, 1] ~ dbern(psi[i])
    for (t in 2:nyear) {
      logit(gamma[i, t-1]) <- alpha.gamma[t-1] + beta.gamma*dist_road[i]
      logit(phi[i, t-1])   <- alpha.phi[t-1]   + beta.phi*dist_road[i]
      muZ[i, t] <- z[i, t-1]*phi[i, t-1] + (1 - z[i, t-1])*gamma[i, t-1]
      z[i, t] ~ dbern(muZ[i, t])
    }
  }

  # Observation process
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

  # Derived
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

model_file <- file.path(OUTPUT_DIR, "model_dynocc_BCAK.txt")
writeLines(model_code, model_file)

# ============================================================================
# PREPARE JAGS INPUT
# ============================================================================

temp_clean    <- d$temp;    temp_clean[is.na(temp_clean)]       <- 0
clutter_clean <- d$clutter; clutter_clean[is.na(clutter_clean)] <- 0
julian_clean  <- d$julian;  julian_clean[is.na(julian_clean)]   <- 0

jags_data <- list(
  y = d$y, J = d$J, nsurv = d$nsurv,
  nsite = d$nsite, nyear = nyear, nregion = d$nregion,
  elevation = d$elevation, water = d$water,
  dist_road = d$dist_road,
  region = d$region,
  temp = temp_clean, clutter = clutter_clean, julian = julian_clean
)

init_z <- apply(d$y, 1:2, function(x) {
  if (all(is.na(x))) return(1)
  ifelse(any(x == 1, na.rm = TRUE), 1, rbinom(1, 1, 0.5))
})

inits <- function() list(
  z = init_z, beta.psi = rnorm(3, 0, 0.5),
  beta.gamma = rnorm(1, 0, 0.5), beta.phi = rnorm(1, 0, 0.5),
  beta.p = rnorm(4, 0, 0.5),
  alpha.gamma = rnorm(nyear - 1, -1, 0.5),
  alpha.phi = rnorm(nyear - 1, 1, 0.5),
  sigma.region = runif(1, 0.1, 1)
)

params <- c("beta.psi", "beta.gamma", "beta.phi", "beta.p",
            "alpha.gamma", "alpha.phi", "gamma.avg", "phi.avg", "ext.avg",
            "sigma.region", "u.region", "psi.fs", "n.occ",
            "mean.gamma", "mean.phi", "mean.ext", "mean.p",
            "lambda", "trend.slope")

# ============================================================================
# FIT
# ============================================================================

cat("Fitting (", N_ITER, " iter, ", N_CHAINS, " chains)...\n", sep = "")
start_time <- Sys.time()

fit <- jags(data = jags_data, inits = inits, parameters.to.save = params,
            model.file = model_file,
            n.chains = N_CHAINS, n.iter = N_ITER, n.burnin = N_BURNIN,
            n.thin = N_THIN, parallel = TRUE)

cat("Done in", round(difftime(Sys.time(), start_time, units = "mins"), 1), "min\n\n")

# ============================================================================
# SAVE & PRINT
# ============================================================================

out_tag <- paste0("QUAD_", SPECIES)
saveRDS(fit, file.path(OUTPUT_DIR, paste0(out_tag, "_fit.rds")))

summ <- fit$summary
summ_df <- as.data.frame(summ); summ_df$parameter <- rownames(summ)
write.csv(summ_df, file.path(OUTPUT_DIR, paste0(out_tag, "_summary.csv")), row.names = FALSE)

sig <- function(lo, hi) ifelse(lo > 0 | hi < 0, " *", "")
pr <- function(label, param) {
  cat(sprintf("  %-20s %7.3f (%7.3f, %7.3f)%s\n",
              label, summ[param, "mean"], summ[param, "2.5%"], summ[param, "97.5%"],
              sig(summ[param, "2.5%"], summ[param, "97.5%"])))
}

cat("RESULTS:", SPECIES, "(BCAK)\n")
cat("Convergence: max Rhat =", round(max(summ[,"Rhat"], na.rm=TRUE), 3),
    "| min n.eff =", round(min(summ[,"n.eff"], na.rm=TRUE)), "\n\n")

cat("Initial Occupancy (psi):\n")
pr("Intercept", "beta.psi[1]"); pr("Elevation", "beta.psi[2]")
pr("Water", "beta.psi[3]"); pr("Region sigma", "sigma.region")

# JAGS names scalar params without index: "beta.gamma" not "beta.gamma[1]"
cat("\nColonization (gamma):\n")
pr("Dist to Road", "beta.gamma")
for (t in 1:(nyear-1))
  pr(paste0(d$years[t], "->", d$years[t+1]), paste0("gamma.avg[", t, "]"))

cat("\nPersistence (phi):\n")
pr("Dist to Road", "beta.phi")
for (t in 1:(nyear-1))
  pr(paste0(d$years[t], "->", d$years[t+1]), paste0("phi.avg[", t, "]"))

cat("\nDetection (p):\n")
pr("Intercept", "beta.p[1]"); pr("Clutter", "beta.p[2]")
pr("Temperature", "beta.p[3]"); pr("Julian", "beta.p[4]")

cat("\nDerived:\n")
pr("Mean colonization", "mean.gamma"); pr("Mean persistence", "mean.phi")
pr("Mean extinction", "mean.ext"); pr("Mean detection", "mean.p")
pr("lambda", "lambda"); pr("Trend slope", "trend.slope")

lam <- summ["lambda", ]
if (lam["97.5%"] < 1) cat("\n  -> DECLINE\n") else
  if (lam["2.5%"] > 1) cat("\n  -> INCREASE\n") else
    cat("\n  -> No clear trend\n")

cat("\nSaved:", file.path(OUTPUT_DIR, paste0(out_tag, "_fit.rds")), "\n")

}
