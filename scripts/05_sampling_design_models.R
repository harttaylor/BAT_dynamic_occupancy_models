################################################################################
# BC Bat Monitoring - Power Analysis (2017-2024) v5
#
# BALANCED GROUP ASSIGNMENT - Equal sites per year within regions
#
# Covariates:
#   ψ (initial): Water Nearby (binary) + Dist to Road
#   γ (colonization): Dist to Road + Dist to Harvest
#   φ (persistence): Dist to Road + Dist to Harvest
#
# REQUIRES: Run 1a_prepare_data_v5.R first
#
################################################################################

library(tidyverse)
library(jagsUI)

# ============================================================================
# SETTINGS
# ============================================================================

SPECIES <- "MYSE"  # Options: LACI, LANO, LABO, MYLU, MYSE
# 01/02/26 launched Laci, lano, labo 

INPUT_DIR <- "data/processed/bc_only_2017start"
OUTPUT_DIR <- paste0("outputs/bc_only_2017start/power_analysis/", SPECIES)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# MCMC settings
n_iter <- 30000
n_burnin <- 15000
n_thin <- 5
n_chains <- 3

# Single iteration per design (set seed for reproducibility)
N_RANDOM_ITERATIONS <- 1
BASE_SEED <- 12345

ROTATING_GRTS <- c(
  186490, 137274, 77162, 243562, 226090, 128058, 114602, 143274,
  16426, 175354, 37274, 283586, 212906, 6202, 172090, 986, 
  265749, 317428, 38442, 132650, 115322, 206394, 26170, 
  230186, 17706, 145578, 58842, 110378
)

years <- 2017:2024
nyear <- length(years)
START_YEAR <- 2017

cat("=" , rep("=", 70), "\n", sep = "")

cat("POWER ANALYSIS (2017-2024) v5:", SPECIES, "\n")
cat("Iterations per design:", N_RANDOM_ITERATIONS, "\n")
cat("=" , rep("=", 70), "\n\n", sep = "")

# ============================================================================
# LOAD DATA
# ============================================================================

model_data <- readRDS(paste0(INPUT_DIR, "/QUAD_", SPECIES, "_model_data.rds"))

y_full <- model_data$y
x.psi <- model_data$x.psi
x.p <- model_data$x.p
x.col <- model_data$x.col
x.per <- model_data$x.per
region <- model_data$region
nregion <- model_data$nregion
nbeta.psi <- model_data$nbeta.psi
nbeta.p <- model_data$nbeta.p
nsite <- model_data$nsite
nvisit <- dim(y_full)[3]
all_sites <- model_data$sites

site_info <- model_data$site_info %>%
  mutate(
    is_rotating = grts_cell %in% ROTATING_GRTS,
    site_idx = row_number(),
    dist_road_std = x.col[, 2],
    dist_harvest_std = x.col[, 3]
  )

n_core <- sum(!site_info$is_rotating)
n_rotating <- sum(site_info$is_rotating)

cat("Sites:", nsite, "(", n_core, "core +", n_rotating, "rotating)\n")
cat("Years:", nyear, "\n")
cat("Regions:", paste(unique(site_info$region), collapse = ", "), "\n\n")

# Show rotating by region
cat("Rotating sites by region:\n")
print(site_info %>% filter(is_rotating) %>% count(region))

# ============================================================================
# BALANCED GROUP ASSIGNMENT
# ============================================================================

assign_balanced_groups <- function(site_info, n_groups, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  assignments <- rep(0, nrow(site_info))
  
  for (reg in unique(site_info$region)) {
    reg_rotating_idx <- which(site_info$is_rotating & site_info$region == reg)
    n_reg <- length(reg_rotating_idx)
    if (n_reg == 0) next
    
    # Random shuffle then round-robin assignment
    shuffled_idx <- sample(reg_rotating_idx)
    for (i in seq_along(shuffled_idx)) {
      assignments[shuffled_idx[i]] <- ((i - 1) %% n_groups) + 1
    }
  }
  return(assignments)
}

# ============================================================================
# MODEL CODE
# ============================================================================

model_code <- "
model {
  # Priors - Initial occupancy
  for(b in 1:nbeta.psi) { beta.psi[b] ~ dnorm(0, 0.368) }
  
  # Region random effect
  sigma.region ~ dunif(0, 5)
  tau.region <- 1 / (sigma.region * sigma.region)
  for(r in 1:nregion) { u.region[r] ~ dnorm(0, tau.region) }
  
  # Colonization
  for(t in 1:(nyear-1)) { alpha.gamma[t] ~ dnorm(0, 0.368) }
  beta.col.road ~ dnorm(0, 0.368)
  beta.col.harvest ~ dnorm(0, 0.368)
  
  # Persistence
  for(t in 1:(nyear-1)) { alpha.phi[t] ~ dnorm(0, 0.368) }
  beta.per.road ~ dnorm(0, 0.368)
  beta.per.harvest ~ dnorm(0, 0.368)
  
  # Detection
  for(b in 1:nbeta.p) { beta.p[b] ~ dnorm(0, 0.368) }
  
  # State process
  for(i in 1:nsite) {
    logit(psi[i]) <- inprod(beta.psi[], x.psi[i,]) + u.region[region[i]]
    z[i, 1] ~ dbern(psi[i])
    
    for(k in 2:nyear) {
      logit(gamma[i, k-1]) <- alpha.gamma[k-1] + beta.col.road * x.col[i, 2] + beta.col.harvest * x.col[i, 3]
      logit(phi[i, k-1]) <- alpha.phi[k-1] + beta.per.road * x.per[i, 2] + beta.per.harvest * x.per[i, 3]
      muZ[i, k] <- z[i, k-1] * phi[i, k-1] + (1 - z[i, k-1]) * gamma[i, k-1]
      z[i, k] ~ dbern(muZ[i, k])
    }
  }
  
  # Observation process
  for(i in 1:nsite) {
    for(k in 1:nyear) {
      for(j in 1:J[i, k]) {
        logit(p[i, k, j]) <- inprod(beta.p[], x.p[i, k, nsurv[i, k, j], ])
        y[i, k, nsurv[i, k, j]] ~ dbern(z[i, k] * p[i, k, j])
      }
    }
  }
  
  # Derived
  for(k in 1:nyear) { psi.fs[k] <- mean(z[, k]) }
  for(t in 1:(nyear-1)) {
    gamma.avg[t] <- ilogit(alpha.gamma[t])
    phi.avg[t] <- ilogit(alpha.phi[t])
  }
  mean.gamma <- mean(gamma.avg[])
  mean.phi <- mean(phi.avg[])
  lambda <- psi.fs[nyear] / max(psi.fs[1], 0.001)
  mean.p <- ilogit(beta.p[1])
}
"

model_file <- paste0(OUTPUT_DIR, "/model.txt")
writeLines(model_code, model_file)

# ============================================================================
# INITIAL VALUES
# ============================================================================

create_inits <- function(y_data) {
  z_init <- apply(y_data, 1:2, function(x) {
    if (all(is.na(x))) return(NA)
    if (any(x == 1, na.rm = TRUE)) return(1)
    return(0)
  })
  
  list(
    z = z_init,
    beta.psi = c(0, 0, 0),
    beta.p = c(0, 0, 0, 0),
    beta.col.road = 0,
    beta.col.harvest = 0,
    beta.per.road = 0,
    beta.per.harvest = 0,
    alpha.gamma = rep(-1, nyear - 1),
    alpha.phi = rep(2, nyear - 1),
    sigma.region = 0.5
  )
}

# ============================================================================
# FIT FUNCTION
# ============================================================================

fit_scenario <- function(scenario_name, n_year_rotation = NULL, groups = NULL, 
                         iteration = 1, save_fit = FALSE) {
  
  cat(sprintf("\n  %s iter %d...", scenario_name, iteration))
  
  if (is.null(n_year_rotation)) {
    y <- y_full
    sampled <- matrix(TRUE, nsite, nyear)
  } else {
    sampled <- matrix(FALSE, nsite, nyear)
    for (i in 1:nsite) {
      for (k in 1:nyear) {
        if (groups[i] == 0) {
          sampled[i, k] <- TRUE
        } else {
          active_group <- ((years[k] - START_YEAR) %% n_year_rotation) + 1
          sampled[i, k] <- (groups[i] == active_group)
        }
      }
    }
    y <- y_full
    for (i in 1:nsite) {
      for (k in 1:nyear) {
        if (!sampled[i, k]) y[i, k, ] <- NA
      }
    }
  }
  
  J <- apply(y, 1:2, function(x) sum(!is.na(x)))
  nsurv <- array(NA, dim = c(nsite, nyear, nvisit))
  for (i in 1:nsite) {
    for (k in 1:nyear) {
      valid <- which(!is.na(y[i, k, ]))
      if (length(valid) > 0) nsurv[i, k, 1:length(valid)] <- valid
    }
  }
  
  n_site_years <- sum(sampled)
  n_detections <- sum(y == 1, na.rm = TRUE)
  sites_per_year <- colSums(sampled)
  
  jags_data <- list(
    y = y, J = J, nsurv = nsurv,
    nsite = nsite, nyear = nyear,
    x.psi = x.psi, nbeta.psi = nbeta.psi,
    x.p = x.p, nbeta.p = nbeta.p,
    x.col = x.col, x.per = x.per,
    region = region, nregion = nregion
  )
  
  inits <- create_inits(y)
  inits_func <- function() {
    inits$sigma.region <- 0.5 + runif(1, -0.01, 0.01)
    return(inits)
  }
  
  params <- c("beta.psi", "beta.p", "beta.col.road", "beta.col.harvest",
              "beta.per.road", "beta.per.harvest", "gamma.avg", "phi.avg", 
              "psi.fs", "mean.gamma", "mean.phi", "mean.p", "lambda", "sigma.region")
  
  fit <- jags(
    data = jags_data, inits = inits_func, parameters.to.save = params,
    model.file = model_file,
    n.chains = n_chains, n.iter = n_iter, n.burnin = n_burnin, n.thin = n_thin,
    parallel = TRUE, verbose = FALSE
  )
  
  summ <- fit$summary
  
  if (save_fit) {
    saveRDS(fit, paste0(OUTPUT_DIR, "/fit_", gsub("-", "", tolower(scenario_name)), ".rds"))
  }
  
  results <- data.frame(
    scenario = scenario_name,
    iteration = iteration,
    site_years = n_site_years,
    pct_of_full = round(n_site_years / (nsite * nyear) * 100, 1),
    detections = n_detections,
    sites_yr_min = min(sites_per_year),
    sites_yr_max = max(sites_per_year),
    lambda_mean = summ["lambda", "mean"],
    lambda_lcl = summ["lambda", "2.5%"],
    lambda_ucl = summ["lambda", "97.5%"],
    ci_width = summ["lambda", "97.5%"] - summ["lambda", "2.5%"],
    beta_psi_water = summ["beta.psi[2]", "mean"],
    beta_psi_road = summ["beta.psi[3]", "mean"],
    beta_col_road = summ["beta.col.road", "mean"],
    beta_col_harvest = summ["beta.col.harvest", "mean"],
    beta_per_road = summ["beta.per.road", "mean"],
    beta_per_harvest = summ["beta.per.harvest", "mean"],
    mean_gamma = summ["mean.gamma", "mean"],
    mean_phi = summ["mean.phi", "mean"],
    mean_p = summ["mean.p", "mean"],
    max_rhat = max(summ[, "Rhat"], na.rm = TRUE)
  )
  
  for (k in 1:nyear) {
    results[[paste0("psi_", years[k])]] <- summ[paste0("psi.fs[", k, "]"), "mean"]
  }
  
  cat(sprintf(" λ=%.2f (%.2f,%.2f)", results$lambda_mean, results$lambda_lcl, results$lambda_ucl))
  
  return(results)
}

# ============================================================================
# FIT ALL SCENARIOS
# ============================================================================

all_results <- list()
idx <- 1

# Full
cat("\n", rep("=", 50), "\n", sep = "")
cat("FULL DESIGN\n")
all_results[[idx]] <- fit_scenario("Full", save_fit = TRUE)
idx <- idx + 1

# Rotation designs
for (n_yr in c(2, 3, 4, 5)) {
  scenario <- paste0(n_yr, "-Year")
  cat("\n", rep("=", 50), "\n", sep = "")
  cat(scenario, "ROTATION\n")
  
  for (iter in 1:N_RANDOM_ITERATIONS) {
    seed <- BASE_SEED + n_yr * 100 + iter
    groups <- assign_balanced_groups(site_info, n_yr, seed = seed)
    
    save_this <- (iter == 1)
    all_results[[idx]] <- fit_scenario(scenario, n_yr, groups, iter, save_this)
    idx <- idx + 1
  }
}

# ============================================================================
# SUMMARIZE
# ============================================================================

results_df <- bind_rows(all_results)

# Set factor order so Full comes first
scenario_order <- c("Full", "2-Year", "3-Year", "4-Year", "5-Year")
results_df$scenario <- factor(results_df$scenario, levels = scenario_order)

summary_df <- results_df %>%
  arrange(scenario) %>%
  mutate(
    ci_width = lambda_ucl - lambda_lcl
  )

full_ci <- summary_df$ci_width[summary_df$scenario == "Full"]
summary_df$precision_loss_pct <- round((summary_df$ci_width / full_ci - 1) * 100, 1)
summary_df$effort_saved_pct <- round(100 - summary_df$pct_of_full, 1)

# ============================================================================
# OUTPUT
# ============================================================================

cat("\n\n", rep("=", 70), "\n", sep = "")
cat("RESULTS:", SPECIES, "\n")
cat(rep("=", 70), "\n\n", sep = "")

print(summary_df %>% 
        select(scenario, site_years, effort_saved_pct, lambda_mean, 
               lambda_lcl, lambda_ucl, ci_width, precision_loss_pct, 
               beta_col_road, beta_col_harvest) %>%
        mutate(across(where(is.numeric), ~round(., 3))))

write_csv(results_df, paste0(OUTPUT_DIR, "/results.csv"))

cat("\nSaved to:", OUTPUT_DIR, "\n")

