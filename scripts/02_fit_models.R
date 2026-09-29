################################################################################
# 02_fit_models.R
#
# Dynamic occupancy model. Section 4 fits ONE species (set SPECIES below);
# section 7 fits every species in turn.
#
#   psi[i]     = a.psi[region] + elevation + dist_water
#   gamma[i,t] = a.gam[region] + year effect + elevation + dist_water
#                              + dist_road + dist_harvest[i,t]
#   phi[i,t]   = a.phi[region] + year effect + elevation + dist_water
#                              + dist_road + dist_harvest[i,t]
#   p[i,t,j]   = a.p + year effect + Alaska offset + clutter + temperature + julian
#
# - Region is a fixed effect (4-5 levels are too few for a variance, and pooling
#   fabricates occupancy in regions where the species doesn't occur). Region and
#   the static covariates sit on colonisation and persistence too, otherwise
#   they wash out after 2017.
# - Year effects are random and unstructured. The year effect on detection
#   (eps.p) absorbs changes in the observation process (classifier, hardware,
#   vetting); a linear p trend would be confounded with the occupancy trend.
# - b.p.ak lets detection differ between Alaska and BC sites (different
#   programs, equipment and auto-ID classifiers). In the BC-only run there are
#   no Alaska sites, so b.p.ak has no data and simply returns its prior: ignore
#   it there.
# - A species is assumed ABSENT from any region where it was never detected
#   (present[r] = 0): its occupancy there is fixed at zero. Without this, the
#   model has no data to rule out "occupied in 2017, then lost", so the first
#   year in those regions is set by the prior and shows a spurious decline.
# - Trends are computed post hoc from posterior draws (slope, P(increase)).
# - Covariates are standardised distances: a NEGATIVE slope means the rate is
#   HIGHER CLOSER to the feature.
#
# Input: one file per species, QUAD_<SPECIES>_data.rds, from 01b
################################################################################

library(jagsUI)

# ============================================================================
# SETTINGS
# ============================================================================

# Two analyses use this script, one per report. Keep ONE block active and put
# a # in front of every line of the other. Use the same block in 01a, 01b, 02.

# --- BC report: every species ---
INPUT_DIR   <- "data/processed/BC"
OUTPUT_DIR  <- "outputs/BC/fits"
det_summary <- read.csv(file.path(INPUT_DIR, "species_detection_summary.csv"))
ALL_SPECIES <- det_summary$species

# --- Alaska report: only species detected in Alaska ---
INPUT_DIR   <- "data/processed/BC_AK"
OUTPUT_DIR  <- "outputs/BC_AK/fits"
det_summary <- read.csv(file.path(INPUT_DIR, "species_detection_summary.csv"))
ALL_SPECIES <- det_summary$species[det_summary$det_Alaska > 0]

ALL_SPECIES                    # species the loop in section 7 will fit

SPECIES <- "COTO"              # species fitted by the step-by-step part (section 4)

MIN_DETECTIONS <- 30           # loop only: below this, logged but not fitted

N_ITER   <- 30000               # 30000 for real runs, more for ANPA
N_BURNIN <- 15000               # 15000 for real runs
N_THIN   <- 10                 # thin hard: z is nsite x nyear per draw

N_CHAINS <- 3

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# 1. LOAD THE SPECIES DATA
# ============================================================================

d <- readRDS(file.path(INPUT_DIR, paste0("QUAD_", SPECIES, "_data.rds")))

# Everything the model needs must be present and numeric. JAGS reports a
# missing element as "not numeric", so check here and fail with a clear message.
needed <- c("y", "J", "nsurv", "nsite", "nyear", "nregion", "region",
            "elevation", "dist_water", "dist_road", "dist_harvest", "clutter",
            "temp", "julian", "surv_ind", "n_surv", "reg_ind", "n_reg")
missing <- needed[!needed %in% names(d) | sapply(needed, function(v) is.null(d[[v]]))]
if (length(missing) > 0)
  stop("Missing from the data file: ", paste(missing, collapse = ", "),
       "\nFile contains: ", paste(names(d), collapse = ", "),
       "\nRe-run 01b so it writes these to ", INPUT_DIR)

not_num <- needed[!sapply(needed, function(v) is.numeric(d[[v]]))]
if (length(not_num) > 0)
  stop("Not numeric (strip units/factor classes in 01b): ",
       paste(not_num, collapse = ", "))

nsite   <- d$nsite
nyear   <- d$nyear
nregion <- d$nregion

# Regions where the species was detected at least once (1) or never (0).
# Occupancy is fixed at zero where present = 0 (see header).
present <- rep(0, nregion)
for (r in 1:nregion) present[r] <- as.numeric(any(d$y[d$region == r, , ] == 1, na.rm = TRUE))
present

# ============================================================================
# 2. JAGS MODEL
# ============================================================================
# Priors: dnorm(0, 0.368) is precision 0.368, i.e. SD 1.65 on the logit scale.
# Approximately flat on the probability scale for an intercept, weakly
# informative for a slope on a standardised covariate, and it stops a
# zero-detection region's intercept running to minus infinity.
# Year-effect SDs get a half-Cauchy(0, 2.5), better than a uniform when the
# true variance is near zero.

model_code <- "
model {

  # ---------------- initial occupancy (2017) ----------------
  for (r in 1:nregion) { a.psi[r] ~ dnorm(0, 0.368) }
  b.psi.elev  ~ dnorm(0, 0.368)
  b.psi.water ~ dnorm(0, 0.368)

  # ---------------- colonisation ----------------
  for (r in 1:nregion) { a.gam[r] ~ dnorm(0, 0.368) }
  b.gam.elev  ~ dnorm(0, 0.368)
  b.gam.water ~ dnorm(0, 0.368)
  b.gam.road  ~ dnorm(0, 0.368)
  b.gam.harv  ~ dnorm(0, 0.368)

  sigma.gam ~ dt(0, 0.16, 1) T(0, )
  tau.gam <- 1 / (sigma.gam * sigma.gam)
  for (t in 1:(nyear - 1)) { eps.gam[t] ~ dnorm(0, tau.gam) }

  # ---------------- persistence ----------------
  for (r in 1:nregion) { a.phi[r] ~ dnorm(0, 0.368) }
  b.phi.elev  ~ dnorm(0, 0.368)
  b.phi.water ~ dnorm(0, 0.368)
  b.phi.road  ~ dnorm(0, 0.368)
  b.phi.harv  ~ dnorm(0, 0.368)

  sigma.phi ~ dt(0, 0.16, 1) T(0, )
  tau.phi <- 1 / (sigma.phi * sigma.phi)
  for (t in 1:(nyear - 1)) { eps.phi[t] ~ dnorm(0, tau.phi) }

  # ---------------- detection ----------------
  a.p         ~ dnorm(0, 0.368)
  b.p.ak      ~ dnorm(0, 0.368)   # Alaska vs BC difference in detection
  b.p.clutter ~ dnorm(0, 0.368)
  b.p.temp    ~ dnorm(0, 0.368)
  b.p.julian  ~ dnorm(0, 0.368)

  sigma.p ~ dt(0, 0.16, 1) T(0, )
  tau.p <- 1 / (sigma.p * sigma.p)
  for (t in 1:nyear) { eps.p[t] ~ dnorm(0, tau.p) }

  # ---------------- ecological process ----------------
  for (i in 1:nsite) {

    logit(psi[i]) <- a.psi[region[i]] +
                     b.psi.elev  * elevation[i] +
                     b.psi.water * dist_water[i]

    z[i, 1] ~ dbern(psi[i] * present[region[i]])

    for (t in 2:nyear) {

      logit(gamma[i, t-1]) <- a.gam[region[i]] + eps.gam[t-1] +
                              b.gam.elev  * elevation[i] +
                              b.gam.water * dist_water[i] +
                              b.gam.road  * dist_road[i] +
                              b.gam.harv  * dist_harvest[i, t-1]

      logit(phi[i, t-1])   <- a.phi[region[i]] + eps.phi[t-1] +
                              b.phi.elev  * elevation[i] +
                              b.phi.water * dist_water[i] +
                              b.phi.road  * dist_road[i] +
                              b.phi.harv  * dist_harvest[i, t-1]

      muZ[i, t] <- present[region[i]] *
                   (z[i, t-1] * phi[i, t-1] + (1 - z[i, t-1]) * gamma[i, t-1])
      z[i, t] ~ dbern(muZ[i, t])
    }
  }

  # ---------------- observation process ----------------
  # J[i,t] = 0 for unsurveyed site-years; those contribute nothing.
  for (i in 1:nsite) {
    for (t in 1:nyear) {
      for (j in 1:J[i, t]) {
        logit(p[i, t, j]) <- a.p + eps.p[t] +
                             b.p.ak      * is_ak[i] +
                             b.p.clutter * clutter[i] +
                             b.p.temp    * temp[i, t, nsurv[i, t, j]] +
                             b.p.julian  * julian[i, t, nsurv[i, t, j]]
        y[i, t, nsurv[i, t, j]] ~ dbern(z[i, t] * p[i, t, j])
      }
    }
  }

  # ---------------- derived: occupancy across all sites ----------------
  # psi.fs  : all sites, including unsurveyed site-years (carried forward by
  #           the Markov process)
  # psi.surv: surveyed site-years only. Effort is uneven across years, so the
  #           two can diverge. Report both.
  for (t in 1:nyear) {
    psi.fs[t]   <- mean(z[, t])
    psi.surv[t] <- inprod(z[, t], surv_ind[, t]) / max(n_surv[t], 1)
    n.occ[t]    <- sum(z[, t])
    p.year[t]   <- ilogit(a.p + eps.p[t])          # BC sites
  }
  mean.p <- ilogit(a.p)

  # ---------------- derived: by region ----------------
  # equilib.reg is the occupancy each region's rates imply. A trajectory that
  # simply walks toward it is the Markov chain relaxing, not a trend.
  for (r in 1:nregion) {
    for (t in 1:nyear) {
      n.occ.reg[r, t] <- inprod(z[, t], reg_ind[, r])
      psi.reg[r, t]   <- n.occ.reg[r, t] / n_reg[r]
    }
    for (t in 1:(nyear - 1)) {
      gamma.reg[r, t] <- ilogit(a.gam[r] + eps.gam[t])
      phi.reg[r, t]   <- ilogit(a.phi[r] + eps.phi[t])
      ext.reg[r, t]   <- 1 - phi.reg[r, t]
    }
    mean.gamma.reg[r] <- mean(gamma.reg[r, ])
    mean.phi.reg[r]   <- mean(phi.reg[r, ])
    mean.ext.reg[r]   <- 1 - mean.phi.reg[r]
    equilib.reg[r]    <- mean.gamma.reg[r] / (mean.gamma.reg[r] + mean.ext.reg[r])
    psi.reg.init[r]   <- ilogit(a.psi[r])
  }
}
"

model_file <- file.path(OUTPUT_DIR, "model_dynocc.txt")
writeLines(model_code, model_file)

# ============================================================================
# 3. JAGS INPUT, INITIAL VALUES, PARAMETERS TO MONITOR
# ============================================================================

jags_data <- list(
  y     = d$y,
  J     = d$J,
  nsurv = d$nsurv,
  
  nsite   = nsite,
  nyear   = nyear,
  nregion = nregion,
  
  region       = d$region,
  elevation    = d$elevation,
  dist_water   = d$dist_water,
  dist_road    = d$dist_road,
  dist_harvest = d$dist_harvest,
  clutter      = d$clutter,
  is_ak        = as.numeric(d$regions[d$region] == "Alaska"),   # 1 = Alaska site
  
  temp   = d$temp,
  julian = d$julian,
  
  surv_ind = d$surv_ind,
  n_surv   = d$n_surv,
  reg_ind  = d$reg_ind,
  n_reg    = d$n_reg,
  present  = present
)

# z must start at 1 wherever the species was detected (or the likelihood is
# impossible), at 0 in regions where it was never detected, and is drawn at
# random elsewhere.
init_z <- apply(d$y, c(1, 2), function(x) {
  if (all(is.na(x))) return(1)
  ifelse(any(x == 1, na.rm = TRUE), 1, rbinom(1, 1, 0.5))
})
init_z[present[d$region] == 0, ] <- 0

inits <- function() list(
  z = init_z,
  a.psi = rnorm(nregion, 0, 0.5),
  a.gam = rnorm(nregion, -1, 0.5),
  a.phi = rnorm(nregion,  1, 0.5),
  b.psi.elev = rnorm(1, 0, 0.3), b.psi.water = rnorm(1, 0, 0.3),
  b.gam.elev = rnorm(1, 0, 0.3), b.gam.water = rnorm(1, 0, 0.3),
  b.gam.road = rnorm(1, 0, 0.3), b.gam.harv  = rnorm(1, 0, 0.3),
  b.phi.elev = rnorm(1, 0, 0.3), b.phi.water = rnorm(1, 0, 0.3),
  b.phi.road = rnorm(1, 0, 0.3), b.phi.harv  = rnorm(1, 0, 0.3),
  a.p = rnorm(1, 0, 0.5),
  b.p.ak = rnorm(1, 0, 0.3),
  b.p.clutter = rnorm(1, 0, 0.3),
  b.p.temp    = rnorm(1, 0, 0.3),
  b.p.julian  = rnorm(1, 0, 0.3),
  sigma.gam = runif(1, 0.1, 1),
  sigma.phi = runif(1, 0.1, 1),
  sigma.p   = runif(1, 0.1, 1)
)

# z is monitored so site-level and regional summaries can be computed later
# without refitting. It is by far the largest object, hence N_THIN = 10.
params <- c(
  "a.psi", "b.psi.elev", "b.psi.water",
  "a.gam", "b.gam.elev", "b.gam.water", "b.gam.road", "b.gam.harv",
  "a.phi", "b.phi.elev", "b.phi.water", "b.phi.road", "b.phi.harv",
  "a.p", "b.p.ak", "b.p.clutter", "b.p.temp", "b.p.julian",
  "sigma.gam", "sigma.phi", "sigma.p",
  "eps.gam", "eps.phi", "eps.p",
  "psi.fs", "psi.surv", "n.occ",
  "psi.reg", "n.occ.reg", "psi.reg.init",
  "gamma.reg", "phi.reg", "ext.reg",
  "mean.gamma.reg", "mean.phi.reg", "mean.ext.reg", "equilib.reg",
  "mean.p", "p.year",
  "z"
)

# ============================================================================
# 4. FIT ONE SPECIES
# ============================================================================

t0 <- Sys.time()

fit <- jags(data = jags_data, inits = inits, parameters.to.save = params,
            model.file = model_file, n.chains = N_CHAINS, n.iter = N_ITER,
            n.burnin = N_BURNIN, n.thin = N_THIN, parallel = TRUE)

runtime_min <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
runtime_min

saveRDS(fit, file.path(OUTPUT_DIR, paste0("QUAD_", SPECIES, "_fit.rds")))

summ <- fit$summary
write.csv(data.frame(parameter = rownames(summ), summ, check.names = FALSE),
          file.path(OUTPUT_DIR, paste0("QUAD_", SPECIES, "_summary.csv")),
          row.names = FALSE)

# ============================================================================
# 5. CONVERGENCE
# ============================================================================
# Judge the parameters that carry the results. z is binary, so its Rhat is not
# meaningful. Intercepts of a zero-detection region are prior-dominated by
# design; judge those on psi.reg instead.

key <- grepl("^(b\\.|a\\.p$|sigma|psi\\.fs|psi\\.surv|psi\\.reg\\[|mean\\.p)", rownames(summ))
not_z <- !grepl("^z\\[", rownames(summ))

convergence <- data.frame(
  max_rhat_key = round(max(summ[key, "Rhat"], na.rm = TRUE), 3),
  min_neff_key = round(min(summ[key, "n.eff"], na.rm = TRUE)),
  max_rhat_all = round(max(summ[not_z, "Rhat"], na.rm = TRUE), 3),
  worst_param  = names(which.max(summ[not_z, "Rhat"])))
convergence

# parameters above Rhat 1.1, if any
round(summ[key & summ[, "Rhat"] > 1.1, c("mean", "2.5%", "97.5%", "Rhat", "n.eff"), drop = FALSE], 3)

# ============================================================================
# 6. RESULTS
# ============================================================================

year_c    <- (1:nyear) - (nyear + 1) / 2
psi_draws <- fit$sims.list$psi.fs
slope_d   <- as.vector(psi_draws %*% year_c / sum(year_c^2))
change_d  <- psi_draws[, nyear] - psi_draws[, 1]

psi_by_year <- data.frame(
  year       = d$years,
  psi_all    = round(summ[paste0("psi.fs[",   1:nyear, "]"), "mean"], 3),
  psi_all_lo = round(summ[paste0("psi.fs[",   1:nyear, "]"), "2.5%"], 3),
  psi_all_hi = round(summ[paste0("psi.fs[",   1:nyear, "]"), "97.5%"], 3),
  psi_surv   = round(summ[paste0("psi.surv[", 1:nyear, "]"), "mean"], 3),
  n_surveyed = d$n_surv)
psi_by_year

trend <- data.frame(
  slope      = round(mean(slope_d), 4),
  slope_lo   = round(unname(quantile(slope_d, 0.025)), 4),
  slope_hi   = round(unname(quantile(slope_d, 0.975)), 4),
  change     = round(mean(change_d), 3),
  p_increase = round(mean(change_d > 0), 3))
trend

# By region. Watch psi_first in a region with zero detections: near zero means
# the fixed region effects are working, above ~0.05 means they are not.
psi_reg_draws <- fit$sims.list$psi.reg      # [draws, region, year]

by_region <- do.call(rbind, lapply(1:nregion, function(r) {
  slp <- as.vector(psi_reg_draws[, r, ] %*% year_c / sum(year_c^2))
  chg <- psi_reg_draws[, r, nyear] - psi_reg_draws[, r, 1]
  data.frame(
    region       = d$regions[r],
    n_sites      = d$n_reg[r],
    n_detections = sum(d$y[d$region == r, , ] == 1, na.rm = TRUE),
    psi_first    = round(mean(psi_reg_draws[, r, 1]), 3),
    psi_last     = round(mean(psi_reg_draws[, r, nyear]), 3),
    slope        = round(mean(slp), 4),
    slope_lo     = round(unname(quantile(slp, 0.025)), 4),
    slope_hi     = round(unname(quantile(slp, 0.975)), 4),
    p_increase   = round(mean(chg > 0), 3),
    colonisation = round(summ[paste0("mean.gamma.reg[", r, "]"), "mean"], 4),
    persistence  = round(summ[paste0("mean.phi.reg[",   r, "]"), "mean"], 3),
    equilibrium  = round(summ[paste0("equilib.reg[",    r, "]"), "mean"], 3))
}))
by_region

# Covariate effects. Negative slopes on the distance covariates mean the rate
# is higher closer to the feature.
round(summ[grepl("^b\\.", rownames(summ)), c("mean", "2.5%", "97.5%", "Rhat")], 3)

# Detection by year (BC sites; Alaska sites add b.p.ak on the logit scale).
# The p trend is a diagnostic, not proof: occupancy and detection trade off.
p_by_year <- data.frame(
  year  = d$years,
  p     = round(summ[paste0("p.year[", 1:nyear, "]"), "mean"], 3),
  p_lo  = round(summ[paste0("p.year[", 1:nyear, "]"), "2.5%"], 3),
  p_hi  = round(summ[paste0("p.year[", 1:nyear, "]"), "97.5%"], 3))
p_by_year

round(c(mean_p = summ["mean.p", "mean"], sigma_p = summ["sigma.p", "mean"],
        p_slope = mean(as.vector(fit$sims.list$p.year %*% year_c / sum(year_c^2)))), 4)

# ============================================================================
# 7. ALL SPECIES
# ============================================================================
# Fits every species in ALL_SPECIES and overwrites any existing fit. Each
# species takes about as long as the single fit in section 4, and nothing
# appears on screen while JAGS runs, so the loop prints when each species
# starts. run_log_species.csv in OUTPUT_DIR is updated after every species,
# so you can open it to see how far the loop has got.

species_log <- NULL
region_log  <- NULL

for (sp in ALL_SPECIES) {
  
  d <- readRDS(file.path(INPUT_DIR, paste0("QUAD_", sp, "_data.rds")))
  nyear <- d$nyear; nregion <- d$nregion
  n_det <- sum(d$y == 1, na.rm = TRUE)
  cat("\n", sp, ": ", n_det, " detections, started ", format(Sys.time(), "%H:%M"), "\n", sep = "")
  
  row <- data.frame(species = sp, n_detections = n_det, n_iter = N_ITER,
                    status = "skipped_too_few", runtime_min = NA)
  
  if (n_det >= MIN_DETECTIONS) {
    
    jags_data$y     <- d$y
    jags_data$J     <- d$J
    jags_data$nsurv <- d$nsurv
    
    present <- rep(0, nregion)
    for (r in 1:nregion) present[r] <- as.numeric(any(d$y[d$region == r, , ] == 1, na.rm = TRUE))
    jags_data$present <- present
    
    init_z <- apply(d$y, c(1, 2), function(x) {
      if (all(is.na(x))) return(1)
      ifelse(any(x == 1, na.rm = TRUE), 1, rbinom(1, 1, 0.5))
    })
    init_z[present[d$region] == 0, ] <- 0
    
    t0  <- Sys.time()
    fit <- tryCatch(
      jags(data = jags_data, inits = inits, parameters.to.save = params,
           model.file = model_file, n.chains = N_CHAINS, n.iter = N_ITER,
           n.burnin = N_BURNIN, n.thin = N_THIN, parallel = TRUE, verbose = FALSE),
      error = function(e) e)
    row$runtime_min <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
    
    if (inherits(fit, "error")) {
      row$status <- paste("failed:", conditionMessage(fit))
    } else {
      row$status <- "ok"
      saveRDS(fit, file.path(OUTPUT_DIR, paste0("QUAD_", sp, "_fit.rds")))
      summ <- fit$summary
      write.csv(data.frame(parameter = rownames(summ), summ, check.names = FALSE),
                file.path(OUTPUT_DIR, paste0("QUAD_", sp, "_summary.csv")),
                row.names = FALSE)
      
      key <- grepl("^(b\\.|a\\.p$|sigma|psi\\.fs|psi\\.surv|psi\\.reg\\[|mean\\.p)", rownames(summ))
      
      year_c   <- (1:nyear) - (nyear + 1) / 2
      psi_d    <- fit$sims.list$psi.fs
      slope_d  <- as.vector(psi_d %*% year_c / sum(year_c^2))
      
      row$max_rhat_key <- round(max(summ[key, "Rhat"], na.rm = TRUE), 3)
      row$min_neff_key <- round(min(summ[key, "n.eff"], na.rm = TRUE))
      row$psi_first    <- round(mean(psi_d[, 1]), 3)
      row$psi_last     <- round(mean(psi_d[, nyear]), 3)
      row$slope        <- round(mean(slope_d), 4)
      row$slope_lo     <- round(unname(quantile(slope_d, 0.025)), 4)
      row$slope_hi     <- round(unname(quantile(slope_d, 0.975)), 4)
      row$p_increase   <- round(mean(psi_d[, nyear] > psi_d[, 1]), 3)
      row$mean_p       <- round(summ["mean.p", "mean"], 3)
      
      psi_reg  <- fit$sims.list$psi.reg
      psi_site <- apply(fit$sims.list$z, c(2, 3), mean)
      ever_det <- apply(d$y, 1, function(x) any(x == 1, na.rm = TRUE))
      
      region_log <- rbind(region_log, do.call(rbind, lapply(1:nregion, function(r) {
        slp   <- as.vector(psi_reg[, r, ] %*% year_c / sum(year_c^2))
        never <- which(d$region == r & !ever_det)
        data.frame(
          species = sp, region = d$regions[r], n_sites = d$n_reg[r],
          n_detections = sum(d$y[d$region == r, , ] == 1, na.rm = TRUE),
          psi_first = round(mean(psi_reg[, r, 1]), 3),
          psi_last  = round(mean(psi_reg[, r, nyear]), 3),
          slope     = round(mean(slp), 4),
          slope_lo  = round(unname(quantile(slp, 0.025)), 4),
          slope_hi  = round(unname(quantile(slp, 0.975)), 4),
          p_increase = round(mean(psi_reg[, r, nyear] > psi_reg[, r, 1]), 3),
          equilibrium = round(summ[paste0("equilib.reg[", r, "]"), "mean"], 3),
          n_never_detected_sites = length(never),
          psi_at_never_detected  = if (length(never)) round(mean(psi_site[never, ]), 3) else NA)
      })))
      
      rm(fit); gc()
    }
  }
  
  # skipped and failed rows have fewer columns, so align before binding
  if (!is.null(species_log)) {
    for (col in setdiff(names(species_log), names(row))) row[[col]] <- NA
    for (col in setdiff(names(row), names(species_log))) species_log[[col]] <- NA
    row <- row[names(species_log)]
  }
  species_log <- rbind(species_log, row)
  
  write.csv(species_log, file.path(OUTPUT_DIR, "run_log_species.csv"), row.names = FALSE)
  if (!is.null(region_log))
    write.csv(region_log, file.path(OUTPUT_DIR, "run_log_regions.csv"), row.names = FALSE)
}

