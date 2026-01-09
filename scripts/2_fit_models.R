################################################################################
# BC Bat Monitoring - Dynamic Occupancy Model Fitting (JAGS)
# 
# COVARIATE STRUCTURE:
#   - Initial occupancy (psi): water [i]
#   - Col/Pers (gamma/phi): clutter + water [i, k-1]
#   - Detection (p): clutter + temp + julian + moon [i, k, nsurv[i,k,j]]
#
################################################################################

library(tidyverse)
library(jagsUI)
library(ggplot2)

# ===== SETTINGS =====

test_species <- "MYLU"

n_iter <- 1000
n_burnin <- 500
n_thin <- 5
n_chains <- 3

# ===== LOAD DATA =====

cat("=== FITTING MODELS FOR", test_species, "===\n\n")

model_data <- readRDS(paste0("data/processed/", test_species, "_model_data.rds"))

# Extract arrays
y <- model_data$y
J <- model_data$J
nsurv <- model_data$nsurv
ind <- model_data$ind
nsite <- model_data$nsite
nyear <- model_data$nyear

x.psi <- model_data$x.psi
x.phi <- model_data$x.phi
x.gamma <- model_data$x.gamma
x.p <- model_data$x.p

nbeta.psi <- model_data$nbeta.psi
nbeta.phi <- model_data$nbeta.phi
nbeta.gamma <- model_data$nbeta.gamma
nbeta.p <- model_data$nbeta.p

years <- model_data$years

cat("Data:", nsite, "sites x", nyear, "years\n")
cat("Observations:", sum(!is.na(y)), "\n")
cat("Detections:", sum(y == 1, na.rm = TRUE), "\n\n")

# ===== MODEL 1: NULL =====

cat("=== MODEL 1: NULL ===\n")

null_model <- "
model{
  # Priors
  psi ~ dunif(0, 1)
  for(t in 1:(nyear-1)){
    phi[t] ~ dunif(0, 1)
    gamma[t] ~ dunif(0, 1)
  }
  p ~ dunif(0, 1)
  
  # State process
  for(i in 1:nsite){
    z[i, 1] ~ dbern(psi)
    for(k in 2:nyear){
      muZ[i, k] <- z[i, k-1]*phi[k-1] + (1 - z[i, k-1])*gamma[k-1]
      z[i, k] ~ dbern(muZ[i, k])
    }
  }
  
  # Observation process
  for(i in 1:nsite){
    for(k in 1:nyear){
      for(j in 1:J[i, k]){
        muy[i, k, j] <- z[i, k]*p
        y[i, k, nsurv[i, k, j]] ~ dbern(muy[i, k, j])
      }
    }
  }
  
  # Derived - use mean(z) for consistency with covariate model
  for(k in 1:nyear){
    psi.fs[k] <- mean(z[, k])
    n.occ[k] <- sum(z[1:nsite, k])
  }
  
  # Extinction rate (for plotting)
  for(t in 1:(nyear-1)){
    ext[t] <- 1 - phi[t]
  }
}
"

writeLines(null_model, "null_model.txt")

jags_null <- list(
  y = y, nsite = nsite, nyear = nyear, J = J, nsurv = nsurv
)

inits_null <- function() {
  list(
    z = apply(y, c(1, 2), function(x) {
      if(all(is.na(x))) return(1)
      return(ifelse(any(x == 1, na.rm = TRUE), 1, rbinom(1, 1, 0.5)))
    }),
    psi = runif(1, 0.3, 0.7),
    phi = runif(nyear - 1, 0.5, 0.9),
    gamma = runif(nyear - 1, 0.05, 0.3),
    p = runif(1, 0.3, 0.7)
  )
}

null_out <- jags(
  data = jags_null,
  inits = inits_null,
  parameters.to.save = c("psi", "phi", "gamma", "ext", "p", "psi.fs", "n.occ"),
  model.file = "null_model.txt",
  n.chains = n_chains,
  n.iter = n_iter, n.burnin = n_burnin, n.thin = n_thin,
  parallel = TRUE
)

cat("\nNULL Model Results:\n")
cat("  Occupancy:", round(null_out$summary["psi", "mean"], 3), "\n")
cat("  Detection:", round(null_out$summary["p", "mean"], 3), "\n")

print(null_out, dig = 3)
saveRDS(null_out, paste0("outputs/", test_species, "_null.rds"))

# ===== MODEL 2: FULL COVARIATES =====

cat("\n=== MODEL 2: FULL COVARIATES ===\n")
cat("  Psi: intercept + water\n")
cat("  Phi/Gamma: intercept + clutter + water\n")
cat("  p: intercept + clutter + temp + julian + moon\n\n")

cov_model <- "
model{
  # ===== PRIORS =====
  
  # Initial occupancy: intercept + water
  for(b in 1:nbeta.psi){
    beta.psi[b] ~ dnorm(0, 0.368)
  }
  
  # Colonization: intercept + clutter + water
  for(b in 1:nbeta.gamma){
    beta.gamma[b] ~ dnorm(0, 0.368)
  }
  
  # Persistence: intercept + clutter + water
  for(b in 1:nbeta.phi){
    beta.phi[b] ~ dnorm(0, 0.368)
  }
  
  # Detection: intercept + clutter + temp + julian + moon
  for(b in 1:nbeta.p){
    beta.p[b] ~ dnorm(0, 0.368)
  }
  
  # ===== STATE PROCESS =====
  
  for(i in 1:nsite){
    # Initial occupancy (year 1)
    logit(psi[i]) <- inprod(beta.psi[], x.psi[i, ])
    z[i, 1] ~ dbern(psi[i])
    
    # Dynamics (years 2 to nyear)
    for(k in 2:nyear){
      # Colonization - uses covariates from year k-1
      logit(gamma[i, k-1]) <- inprod(beta.gamma[], x.gamma[i, k-1, ])
      
      # Persistence - uses covariates from year k-1
      logit(phi[i, k-1]) <- inprod(beta.phi[], x.phi[i, k-1, ])
      
      muZ[i, k] <- z[i, k-1]*phi[i, k-1] + (1 - z[i, k-1])*gamma[i, k-1]
      z[i, k] ~ dbern(muZ[i, k])
    }
  }
  
  # ===== OBSERVATION PROCESS =====
  
  for(i in 1:nsite){
    for(k in 1:nyear){
      for(j in 1:J[i, k]){
        logit(p[i, k, j]) <- inprod(beta.p[], x.p[i, k, nsurv[i, k, j], ])
        muy[i, k, j] <- z[i, k] * p[i, k, j]
        y[i, k, nsurv[i, k, j]] ~ dbern(muy[i, k, j])
      }
    }
  }
  
  # ===== DERIVED QUANTITIES =====
  
  # Annual occupancy
  for(k in 1:nyear){
    psi.fs[k] <- mean(z[, k])
    n.occ[k] <- sum(z[, k])
  }
  
  # Mean colonization and extinction by year
  for(k in 1:(nyear-1)){
    mean.gamma[k] <- mean(gamma[, k])
    mean.phi[k] <- mean(phi[, k])
    mean.ext[k] <- 1 - mean.phi[k]
  }
  
  # Mean detection (at average covariate values = 0 for standardized)
  mean.p <- ilogit(beta.p[1])
}
"

writeLines(cov_model, "cov_model.txt")

jags_cov <- list(
  y = y, nsite = nsite, nyear = nyear, J = J, nsurv = nsurv,
  x.psi = x.psi, nbeta.psi = nbeta.psi,
  x.phi = x.phi, nbeta.phi = nbeta.phi,
  x.gamma = x.gamma, nbeta.gamma = nbeta.gamma,
  x.p = x.p, nbeta.p = nbeta.p
)

inits_cov <- function() {
  list(
    z = apply(y, c(1, 2), function(x) {
      if(all(is.na(x))) return(1)
      return(ifelse(any(x == 1, na.rm = TRUE), 1, rbinom(1, 1, 0.5)))
    })
  )
}

params_cov <- c("beta.psi", "beta.gamma", "beta.phi", "beta.p",
                "psi.fs", "n.occ", "mean.gamma", "mean.phi", "mean.ext", "mean.p",
                "z", "p")

cov_out <- jags(
  data = jags_cov,
  inits = inits_cov,
  parameters.to.save = params_cov,
  model.file = "cov_model.txt",
  n.chains = n_chains,
  n.iter = n_iter, n.burnin = n_burnin, n.thin = n_thin,
  parallel = TRUE
)

print(cov_out, dig = 3)
saveRDS(cov_out, paste0("outputs/", test_species, "_cov.rds"))


# ===== PARAMETER ESTIMATES =====

cat("\n=== PARAMETER ESTIMATES ===\n")

sig_star <- function(lcl, ucl) ifelse(lcl > 0 | ucl < 0, "*", "")

# Initial occupancy
cat("\nINITIAL OCCUPANCY (psi):\n")
cat("  Intercept:  ", sprintf("%6.3f (%6.3f, %6.3f)", 
                              cov_out$summary["beta.psi[1]", "mean"],
                              cov_out$summary["beta.psi[1]", "2.5%"],
                              cov_out$summary["beta.psi[1]", "97.5%"]), "\n")
cat("  Water:      ", sprintf("%6.3f (%6.3f, %6.3f) %s", 
                              cov_out$summary["beta.psi[2]", "mean"],
                              cov_out$summary["beta.psi[2]", "2.5%"],
                              cov_out$summary["beta.psi[2]", "97.5%"],
                              sig_star(cov_out$summary["beta.psi[2]", "2.5%"], cov_out$summary["beta.psi[2]", "97.5%"])), "\n")

# Colonization
cat("\nCOLONIZATION (gamma):\n")
cat("  Intercept:  ", sprintf("%6.3f (%6.3f, %6.3f)", 
                              cov_out$summary["beta.gamma[1]", "mean"],
                              cov_out$summary["beta.gamma[1]", "2.5%"],
                              cov_out$summary["beta.gamma[1]", "97.5%"]), "\n")
cat("  Clutter:    ", sprintf("%6.3f (%6.3f, %6.3f) %s", 
                              cov_out$summary["beta.gamma[2]", "mean"],
                              cov_out$summary["beta.gamma[2]", "2.5%"],
                              cov_out$summary["beta.gamma[2]", "97.5%"],
                              sig_star(cov_out$summary["beta.gamma[2]", "2.5%"], cov_out$summary["beta.gamma[2]", "97.5%"])), "\n")
cat("  Water:      ", sprintf("%6.3f (%6.3f, %6.3f) %s", 
                              cov_out$summary["beta.gamma[3]", "mean"],
                              cov_out$summary["beta.gamma[3]", "2.5%"],
                              cov_out$summary["beta.gamma[3]", "97.5%"],
                              sig_star(cov_out$summary["beta.gamma[3]", "2.5%"], cov_out$summary["beta.gamma[3]", "97.5%"])), "\n")

# Persistence
cat("\nPERSISTENCE (phi):\n")
cat("  Intercept:  ", sprintf("%6.3f (%6.3f, %6.3f)", 
                              cov_out$summary["beta.phi[1]", "mean"],
                              cov_out$summary["beta.phi[1]", "2.5%"],
                              cov_out$summary["beta.phi[1]", "97.5%"]), "\n")
cat("  Clutter:    ", sprintf("%6.3f (%6.3f, %6.3f) %s", 
                              cov_out$summary["beta.phi[2]", "mean"],
                              cov_out$summary["beta.phi[2]", "2.5%"],
                              cov_out$summary["beta.phi[2]", "97.5%"],
                              sig_star(cov_out$summary["beta.phi[2]", "2.5%"], cov_out$summary["beta.phi[2]", "97.5%"])), "\n")
cat("  Water:      ", sprintf("%6.3f (%6.3f, %6.3f) %s", 
                              cov_out$summary["beta.phi[3]", "mean"],
                              cov_out$summary["beta.phi[3]", "2.5%"],
                              cov_out$summary["beta.phi[3]", "97.5%"],
                              sig_star(cov_out$summary["beta.phi[3]", "2.5%"], cov_out$summary["beta.phi[3]", "97.5%"])), "\n")

# Detection
cat("\nDETECTION (p):\n")
cat("  Intercept:  ", sprintf("%6.3f", cov_out$summary["beta.p[1]", "mean"]), "\n")
cat("  Clutter:    ", sprintf("%6.3f (%6.3f, %6.3f) %s", 
                              cov_out$summary["beta.p[2]", "mean"],
                              cov_out$summary["beta.p[2]", "2.5%"],
                              cov_out$summary["beta.p[2]", "97.5%"],
                              sig_star(cov_out$summary["beta.p[2]", "2.5%"], cov_out$summary["beta.p[2]", "97.5%"])), "\n")
cat("  Temp:       ", sprintf("%6.3f (%6.3f, %6.3f) %s", 
                              cov_out$summary["beta.p[3]", "mean"],
                              cov_out$summary["beta.p[3]", "2.5%"],
                              cov_out$summary["beta.p[3]", "97.5%"],
                              sig_star(cov_out$summary["beta.p[3]", "2.5%"], cov_out$summary["beta.p[3]", "97.5%"])), "\n")
cat("  Julian:     ", sprintf("%6.3f (%6.3f, %6.3f) %s", 
                              cov_out$summary["beta.p[4]", "mean"],
                              cov_out$summary["beta.p[4]", "2.5%"],
                              cov_out$summary["beta.p[4]", "97.5%"],
                              sig_star(cov_out$summary["beta.p[4]", "2.5%"], cov_out$summary["beta.p[4]", "97.5%"])), "\n")
cat("  Moon:       ", sprintf("%6.3f (%6.3f, %6.3f) %s", 
                              cov_out$summary["beta.p[5]", "mean"],
                              cov_out$summary["beta.p[5]", "2.5%"],
                              cov_out$summary["beta.p[5]", "97.5%"],
                              sig_star(cov_out$summary["beta.p[5]", "2.5%"], cov_out$summary["beta.p[5]", "97.5%"])), "\n")
cat("  Mean p:     ", sprintf("%6.3f", cov_out$summary["mean.p", "mean"]), "\n")

cat("\n  * = 95% CI excludes 0\n")

# ===== BAYESIAN P-VALUE =====

calculate_bpv <- function(jags_out, y, J, nsurv, n_sims = 500) {
  
  z_samples <- jags_out$sims.list$z
  p_samples <- jags_out$sims.list$p
  
  n_iter <- dim(z_samples)[1]
  nsite <- dim(z_samples)[2]
  nyear <- dim(z_samples)[3]
  
  use_iter <- sample(1:n_iter, min(n_sims, n_iter))
  
  fit_obs <- numeric(length(use_iter))
  fit_sim <- numeric(length(use_iter))
  
  for(s in seq_along(use_iter)) {
    iter <- use_iter[s]
    T_obs <- 0
    T_sim <- 0
    
    for(i in 1:nsite) {
      for(k in 1:nyear) {
        if(J[i, k] > 0) {
          for(j in 1:J[i, k]) {
            idx <- nsurv[i, k, j]
            e_ijk <- z_samples[iter, i, k] * p_samples[iter, i, k, j]
            y_obs <- y[i, k, idx]
            
            if(!is.na(y_obs)) {
              T_obs <- T_obs + (sqrt(y_obs) - sqrt(e_ijk))^2
              y_sim <- rbinom(1, 1, e_ijk)
              T_sim <- T_sim + (sqrt(y_sim) - sqrt(e_ijk))^2
            }
          }
        }
      }
    }
    fit_obs[s] <- T_obs
    fit_sim[s] <- T_sim
  }
  
  bpv <- mean(fit_sim > fit_obs)
  return(list(bpv = bpv, fit_obs = fit_obs, fit_sim = fit_sim))
}

cat("\nCalculating Bayesian p-value...\n")
bpv_result <- calculate_bpv(cov_out, y, J, nsurv)
cat("Bayesian p-value:", round(bpv_result$bpv, 3), "\n")

if(bpv_result$bpv > 0.1 & bpv_result$bpv < 0.9) {
  cat("  Interpretation: Adequate model fit\n")
} else {
  cat("  Interpretation: Potential lack of fit - investigate further\n")
}

saveRDS(bpv_result, paste0("output/", test_species, "_bpv.rds"))

# ===== PLOT 1: OCCUPANCY TREND (from null model) =====

palette.colors(palette = "Dark 2")
palette.pals()

psi_idx <- grep("^psi.fs\\[", rownames(null_out$summary))

trend_df <- data.frame(
  year = years,
  psi = null_out$summary[psi_idx, "mean"],
  lcl = null_out$summary[psi_idx, "2.5%"],
  ucl = null_out$summary[psi_idx, "97.5%"]
)

p1 <- ggplot(trend_df, aes(x = year, y = psi)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.3, fill = "#1B9E77") +
  geom_line(linewidth = 1.2, color = "#1B9E77") +
  geom_point(size = 3, color = "#1B9E77") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = years) +
  labs(x = "Year", y = "Occupancy probability",
       title = paste(test_species, "- Occupancy Trend")) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

print(p1)
ggsave(paste0("output/", test_species, "_occupancy_trend.png"), p1, width = 10, height = 6)

# ===== PLOT 2: DYNAMICS (from null model) - PERSISTENCE VERSION =====

# Create year labels for transitions (nyear - 1 transitions)
n_transitions <- nyear - 1
transition_years <- years[1:n_transitions]

# Extract gamma and phi from null model (these are vectors of length nyear-1)
dynamics_df <- data.frame(
  year = transition_years,
  colonization = null_out$summary[paste0("gamma[", 1:n_transitions, "]"), "mean"],
  col_lcl = null_out$summary[paste0("gamma[", 1:n_transitions, "]"), "2.5%"],
  col_ucl = null_out$summary[paste0("gamma[", 1:n_transitions, "]"), "97.5%"],
  persistence = null_out$summary[paste0("phi[", 1:n_transitions, "]"), "mean"],
  per_lcl = null_out$summary[paste0("phi[", 1:n_transitions, "]"), "2.5%"],
  per_ucl = null_out$summary[paste0("phi[", 1:n_transitions, "]"), "97.5%"]
)

dynamics_long <- dynamics_df %>%
  pivot_longer(cols = c(colonization, persistence), names_to = "rate", values_to = "value") %>%
  mutate(
    lcl = ifelse(rate == "colonization", col_lcl, per_lcl),
    ucl = ifelse(rate == "colonization", col_ucl, per_ucl),
    rate = factor(rate, levels = c("colonization", "persistence"),
                  labels = c("Colonization", "Persistence"))
  )

p2 <- ggplot(dynamics_long, aes(x = year, y = value, color = rate, fill = rate)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.3, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Colonization" = "#D95F02", "Persistence" = "#7570B3")) +
  scale_fill_manual(values = c("Colonization" = "#D95F02", "Persistence" = "#7570B3")) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = years) +
  labs(x = "Year", y = "Probability", color = NULL, fill = NULL,
       title = paste(test_species, "- Colonization and Persistence Dynamics (Null Model)")) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

print(p2)
ggsave(paste0("output/", test_species, "_dynamics.png"), p2, width = 10, height = 6)


# ===== PLOT 3: COVARIATE EFFECTS (from covariate model) =====

effects_df <- data.frame(
  parameter = c("beta.psi[2]", "beta.gamma[2]", "beta.gamma[3]", 
                "beta.phi[2]", "beta.phi[3]",
                "beta.p[2]", "beta.p[3]", "beta.p[4]", "beta.p[5]"),
  process = c("Initial Occ", "Colonization", "Colonization",
              "Persistence", "Persistence",
              "Detection", "Detection", "Detection", "Detection"),
  covariate = c("Water", "Clutter", "Water", "Clutter", "Water",
                "Clutter", "Temp", "Julian", "Moon")
)

effects_df$mean <- cov_out$summary[effects_df$parameter, "mean"]
effects_df$lcl <- cov_out$summary[effects_df$parameter, "2.5%"]
effects_df$ucl <- cov_out$summary[effects_df$parameter, "97.5%"]
effects_df$significant <- effects_df$lcl > 0 | effects_df$ucl < 0
effects_df$label <- paste(effects_df$process, "-", effects_df$covariate)

p3 <- ggplot(effects_df, aes(x = reorder(label, mean), y = mean, color = significant)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_pointrange(aes(ymin = lcl, ymax = ucl), size = 0.8, linewidth = 1) +
  scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "darkred"),
                     name = "95% CI excludes 0") +
  coord_flip() +
  labs(x = NULL, y = "Effect size (logit scale)",
       title = paste(test_species, "- Covariate Effects")) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

print(p3)
ggsave(paste0("output/", test_species, "_covariate_effects.png"), p3, width = 8, height = 6)

# ===== PLOT 5: REGIONAL TRENDS (if available) =====

if(!is.null(site_info) && "region" %in% names(site_info) && length(unique(site_info$region)) > 1) {
  
  z_mean <- cov_out$mean$z
  
  regional_df <- data.frame()
  
  for(k in 1:nyear) {
    region_occ <- data.frame(
      site = site_info$site,
      region = site_info$region,
      z = z_mean[, k],
      year = years[k]
    ) %>%
      group_by(region, year) %>%
      summarize(
        mean_occ = mean(z),
        n_sites = n(),
        .groups = "drop"
      )
    
    regional_df <- bind_rows(regional_df, region_occ)
  }
  
  p5 <- ggplot(regional_df, aes(x = year, y = mean_occ, color = region)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = years) +
    labs(x = "Year", y = "Mean occupancy probability",
         title = paste(test_species, "- Regional Occupancy Trends"),
         color = "Region") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  print(p5)
  ggsave(paste0("output/", test_species, "_regional_trends.png"), p5, width = 10, height = 6)
  cat("Saved: regional_trends.png\n")
}

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Saved to output/:\n")
cat("  ", test_species, "_null.rds\n")
cat("  ", test_species, "_cov.rds\n")
cat("  ", test_species, "_bpv.rds\n")
cat("  ", test_species, "_occupancy_trend.png\n")
cat("  ", test_species, "_dynamics.png\n")
cat("  ", test_species, "_covariate_effects.png\n")
cat("  ", test_species, "_bpv_plot.png\n")