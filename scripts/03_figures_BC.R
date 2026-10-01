################################################################################
# 03_figures_BC.R
#
# Figures and tables from the BC-only run of 01a, 01b and 02.
#
# PART 1 loops over species and makes the per-species figures and tables.

# PARTS 2-4 combine the tables, make the cross-species figures and write a
#   data dictionary.
#
#
# TABLES (outputs/BC/tables/per_species/, combined as ALL_SPECIES_*.csv)
#   covariate_effects_<SP>.csv
#   occupancy_trajectory_<SP>.csv   BC-wide (3 versions) and regional, by year
#   dynamics_<SP>.csv               colonisation and persistence, BC and regional
#   trends_<SP>.csv                 slope, change and P(increase), BC and regional
#   site_year_occupancy_<SP>.csv    one row per site-year (for maps)
#   site_summary_<SP>.csv           one row per site (for maps)
#
################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)

# ============================================================================
# SETTINGS
# ============================================================================

DATA_DIR <- "data/processed/BC"
FIT_DIR  <- "outputs/BC/fits"
OUT_DIR  <- "outputs/BC"

det_summary     <- read.csv(file.path(DATA_DIR, "species_detection_summary.csv"))
SPECIES_TO_PLOT <- det_summary$species     # species without a fit are skipped
# SPECIES_TO_PLOT <- "MYLU"                # uncomment to test one species

TAB_DIR     <- file.path(OUT_DIR, "tables")
per_sp_dir  <- file.path(OUT_DIR, "tables", "per_species")
all_fig_dir <- file.path(OUT_DIR, "figures", "_all_species")
dir.create(per_sp_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(all_fig_dir, recursive = TRUE, showWarnings = FALSE)

SPECIES_NAMES <- c(
  ANPA = "Pallid bat",               COTO = "Townsend's big-eared bat",
  EPFU = "Big brown bat",            EUMA = "Spotted bat",
  LABO = "Eastern red bat",          LACI = "Hoary bat",
  LANO = "Silver-haired bat",        MYCA = "California myotis",
  MYCI = "Western small-footed myotis", MYEV = "Western long-eared myotis",
  MYLU = "Little brown myotis",      MYSE = "Northern myotis",
  MYTH = "Fringed myotis",           MYVO = "Long-legged myotis",
  MYYU = "Yuma myotis",              PAHE = "Canyon bat",
  TABR = "Brazilian free-tailed bat"
)

LINE_COL <- "#1b9e77"

theme_occ <- theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey95", colour = "grey40"),
        strip.text = element_text(size = 10),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey30"))

# ============================================================================
# PART 1. PER-SPECIES FIGURES AND TABLES
# ============================================================================

for (SPECIES in SPECIES_TO_PLOT) {
  
  # ============================================================================
  # 1. LOAD
  # ============================================================================
  
  fit_file <- file.path(FIT_DIR, paste0("QUAD_", SPECIES, "_fit.rds"))
  if (!file.exists(fit_file)) {
    cat("\n", SPECIES, ": no fit (under 30 detections, or not run yet), skipped\n", sep = "")
    next
  }
  cat("\n----", SPECIES, "----\n")
  
  fit <- readRDS(fit_file)
  d   <- readRDS(file.path(DATA_DIR, paste0("QUAD_", SPECIES, "_data.rds")))
  sl  <- fit$sims.list
  
  sp_title <- ifelse(SPECIES %in% names(SPECIES_NAMES),
                     paste0(SPECIES, " \u2014 ", SPECIES_NAMES[SPECIES]), SPECIES)
  sp_fig_dir <- file.path(OUT_DIR, "figures", SPECIES)
  dir.create(sp_fig_dir, recursive = TRUE, showWarnings = FALSE)
  
  nsite     <- d$nsite
  nyear     <- d$nyear
  nregion   <- d$nregion
  years     <- d$years
  regions   <- d$regions
  site_info <- d$site_info
  n_draws   <- length(sl$mean.p)
  
  # ============================================================================
  # 2. OBSERVED DATA SUMMARIES
  # ============================================================================
  
  surveyed <- d$J > 0                                                   # [site, year]
  detected <- apply(d$y, c(1, 2), function(x) any(x == 1, na.rm = TRUE)) # [site, year]
  detected[!surveyed] <- FALSE
  ever_detected <- apply(detected, 1, any)
  
  n_det_total <- sum(d$y == 1, na.rm = TRUE)
  n_det_reg <- rep(NA, nregion)
  for (r in 1:nregion) n_det_reg[r] <- sum(d$y[d$region == r, , ] == 1, na.rm = TRUE)
  
  # Sites in each region where the species was detected at least once
  n_sites_det_reg <- rep(NA, nregion)
  for (r in 1:nregion) n_sites_det_reg[r] <- sum(ever_detected[d$region == r])
  
  # Panel labels. A "detection night" is one survey night (visit) at one site
  # with at least one call identified to this species.
  region_label <- paste0(regions, "\n", d$n_reg, " sites, detected at ", n_sites_det_reg,
                         "\n", n_det_reg, " detection nights")
  region_label[n_det_reg == 0] <- paste0(regions[n_det_reg == 0], "\n",
                                         d$n_reg[n_det_reg == 0], " sites\nNOT DETECTED")
  
  # Naive occupancy: proportion of surveyed sites with at least one detection
  naive_bc <- colSums(detected) / colSums(surveyed)
  
  naive_reg <- matrix(NA, nregion, nyear)
  nsurv_reg <- matrix(NA, nregion, nyear)
  for (r in 1:nregion) {
    rows <- d$region == r
    nsurv_reg[r, ] <- colSums(surveyed[rows, , drop = FALSE])
    naive_reg[r, ] <- ifelse(nsurv_reg[r, ] > 0,
                             colSums(detected[rows, , drop = FALSE]) / nsurv_reg[r, ],
                             NA)
  }
  
  # ============================================================================
  # 3. COVARIATE EFFECTS
  # ============================================================================
  # These slopes are shared across regions, so they ARE the BC-wide effects.
  # All covariates are standardised; distances were log(x+1) first. A slope is
  # the change in logit(probability) per 1 SD of the covariate.
  
  cov_eff <- data.frame(
    param = c("b.psi.elev", "b.psi.water",
              "b.gam.elev", "b.gam.water", "b.gam.road", "b.gam.harv",
              "b.phi.elev", "b.phi.water", "b.phi.road", "b.phi.harv",
              "b.p.clutter", "b.p.temp", "b.p.julian"),
    process = c(rep("Initial occupancy", 2), rep("Colonization", 4),
                rep("Persistence", 4), rep("Detection", 3)),
    covariate = c("Elevation", "Distance to water",
                  "Elevation", "Distance to water", "Distance to road", "Distance to harvest",
                  "Elevation", "Distance to water", "Distance to road", "Distance to harvest",
                  "Clutter", "Temperature", "Julian date"),
    stringsAsFactors = FALSE
  )
  
  cov_eff$mean <- NA; cov_eff$lo95 <- NA; cov_eff$hi95 <- NA
  cov_eff$lo50 <- NA; cov_eff$hi50 <- NA; cov_eff$f <- NA
  
  for (k in 1:nrow(cov_eff)) {
    draws <- sl[[cov_eff$param[k]]]
    cov_eff$mean[k] <- mean(draws)
    cov_eff$lo95[k] <- quantile(draws, 0.025)
    cov_eff$hi95[k] <- quantile(draws, 0.975)
    cov_eff$lo50[k] <- quantile(draws, 0.25)
    cov_eff$hi50[k] <- quantile(draws, 0.75)
    # f = posterior probability the effect has the same sign as its mean
    cov_eff$f[k] <- ifelse(mean(draws) > 0, mean(draws > 0), mean(draws < 0))
  }
  
  cov_eff$excludes_zero <- cov_eff$lo95 > 0 | cov_eff$hi95 < 0
  cov_eff$direction <- ifelse(cov_eff$lo95 > 0, "Positive",
                              ifelse(cov_eff$hi95 < 0, "Negative", "Includes zero"))
  cov_eff$direction <- factor(cov_eff$direction, levels = c("Positive", "Negative", "Includes zero"))
  cov_eff$process <- factor(cov_eff$process,
                            levels = c("Initial occupancy", "Colonization",
                                       "Persistence", "Detection"))
  cov_eff$covariate <- factor(cov_eff$covariate,
                              levels = rev(c("Elevation", "Distance to water",
                                             "Distance to road", "Distance to harvest",
                                             "Clutter", "Temperature", "Julian date")))
  
  fig1 <- ggplot(cov_eff, aes(x = mean, y = covariate, colour = direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_linerange(aes(xmin = lo95, xmax = hi95), linewidth = 1) +
    geom_point(size = 2.5) +
    facet_grid(process ~ ., scales = "free_y", space = "free_y") +
    scale_colour_manual(values = c(Positive = "#2166ac", Negative = "#b2182b",
                                   `Includes zero` = "grey65"),
                        labels = c(Positive = "Positive (95% CI above 0)",
                                   Negative = "Negative (95% CI below 0)",
                                   `Includes zero` = "95% CI includes 0"),
                        name = NULL, drop = FALSE) +
    labs(title = paste0(sp_title, ": covariate effects"),
         subtitle = "Posterior mean and 95% credible interval. BC-wide slopes, shared across regions.",
         x = "Effect on logit scale (per 1 SD of covariate)", y = NULL) +
    theme_occ +
    theme(legend.position = "bottom", strip.text.y = element_text(angle = 0))
  
  ggsave(file.path(sp_fig_dir, paste0("fig1_covariate_effects_", SPECIES, ".png")),
         fig1, width = 8, height = 7, dpi = 300)
  
  # ============================================================================
  # 4. BC-WIDE TRAJECTORIES
  # ============================================================================
  
  z_d <- sl$z                                  # [draws, site, year]
  
  # --- realised colonisation and persistence, all sites ---
  col_bc <- matrix(NA, n_draws, nyear - 1)
  per_bc <- matrix(NA, n_draws, nyear - 1)
  for (t in 2:nyear) {
    prev <- z_d[, , t - 1]
    cur  <- z_d[, , t]
    col_bc[, t - 1] <- rowSums((1 - prev) * cur) / rowSums(1 - prev)
    per_bc[, t - 1] <- rowSums(prev * cur) / rowSums(prev)
  }
  
  col_bc_mean <- colMeans(col_bc, na.rm = TRUE)
  col_bc_lo   <- apply(col_bc, 2, quantile, 0.025, na.rm = TRUE)
  col_bc_hi   <- apply(col_bc, 2, quantile, 0.975, na.rm = TRUE)
  per_bc_mean <- colMeans(per_bc, na.rm = TRUE)
  per_bc_lo   <- apply(per_bc, 2, quantile, 0.025, na.rm = TRUE)
  per_bc_hi   <- apply(per_bc, 2, quantile, 0.975, na.rm = TRUE)
  
  # --- region-weighted BC occupancy: each region counts equally ---
  psi_regw <- apply(sl$psi.reg, c(1, 3), mean)  # [draws, year]
  
  dyn_bc <- rbind(
    data.frame(process = "Occupancy", year = years,
               mean = colMeans(sl$psi.fs),
               lo = apply(sl$psi.fs, 2, quantile, 0.025),
               hi = apply(sl$psi.fs, 2, quantile, 0.975),
               naive = naive_bc),
    data.frame(process = "Colonization", year = years[-1],
               mean = col_bc_mean, lo = col_bc_lo, hi = col_bc_hi, naive = NA),
    data.frame(process = "Persistence", year = years[-1],
               mean = per_bc_mean, lo = per_bc_lo, hi = per_bc_hi, naive = NA)
  )
  dyn_bc$process <- factor(dyn_bc$process, levels = c("Occupancy", "Colonization", "Persistence"))
  
  fig2 <- ggplot(dyn_bc, aes(x = year, y = mean)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = LINE_COL, alpha = 0.2) +
    geom_line(colour = LINE_COL, linewidth = 1) +
    geom_point(colour = LINE_COL, size = 2) +
    facet_wrap(~ process, nrow = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = years) +
    labs(title = paste0(sp_title, ": BC-wide dynamics"),
         subtitle = paste0("Posterior mean (95% CI) across all ", nsite, " sites."),
         x = "Year", y = "Probability") +
    theme_occ +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(sp_fig_dir, paste0("fig2_bc_dynamics_", SPECIES, ".png")),
         fig2, width = 11, height = 4.5, dpi = 300)
  
  # ============================================================================
  # 5. REGIONAL TRAJECTORIES
  # ============================================================================
  
  occ_reg <- data.frame()
  dyn_reg <- data.frame()
  
  for (r in 1:nregion) {
    
    rows <- which(d$region == r)
    psi_r <- sl$psi.reg[, r, ]                   # [draws, year]
    
    occ_reg <- rbind(occ_reg, data.frame(
      region = regions[r], region_label = region_label[r], year = years,
      mean = colMeans(psi_r),
      lo = apply(psi_r, 2, quantile, 0.025),
      hi = apply(psi_r, 2, quantile, 0.975),
      naive = naive_reg[r, ],
      n_surveyed = nsurv_reg[r, ],
      stringsAsFactors = FALSE))
    
    # realised rates within this region
    col_r <- matrix(NA, n_draws, nyear - 1)
    per_r <- matrix(NA, n_draws, nyear - 1)
    for (t in 2:nyear) {
      prev <- z_d[, rows, t - 1, drop = FALSE]
      cur  <- z_d[, rows, t,     drop = FALSE]
      col_r[, t - 1] <- rowSums((1 - prev) * cur) / rowSums(1 - prev)
      per_r[, t - 1] <- rowSums(prev * cur) / rowSums(prev)
    }
    
    col_mean <- colMeans(col_r, na.rm = TRUE)
    col_lo   <- apply(col_r, 2, quantile, 0.025, na.rm = TRUE)
    col_hi   <- apply(col_r, 2, quantile, 0.975, na.rm = TRUE)
    per_mean <- colMeans(per_r, na.rm = TRUE)
    per_lo   <- apply(per_r, 2, quantile, 0.025, na.rm = TRUE)
    per_hi   <- apply(per_r, 2, quantile, 0.975, na.rm = TRUE)
    
    # Blank out rates that are undefined in most draws (e.g. persistence in a
    # region where almost no sites are ever occupied)
    col_ok <- colMeans(is.finite(col_r)) >= 0.5
    per_ok <- colMeans(is.finite(per_r)) >= 0.5
    col_mean[!col_ok] <- NA; col_lo[!col_ok] <- NA; col_hi[!col_ok] <- NA
    per_mean[!per_ok] <- NA; per_lo[!per_ok] <- NA; per_hi[!per_ok] <- NA
    
    # Persistence is also left blank where the species was never detected in
    # the region. No site there is known to have been occupied, so the data
    # contain no information about persistence and any value is just the prior.
    if (n_det_reg[r] == 0) { per_mean[] <- NA; per_lo[] <- NA; per_hi[] <- NA }
    
    dyn_reg <- rbind(dyn_reg,
                     data.frame(region = regions[r], region_label = region_label[r],
                                process = "Occupancy", year = years,
                                mean = colMeans(psi_r),
                                lo = apply(psi_r, 2, quantile, 0.025),
                                hi = apply(psi_r, 2, quantile, 0.975),
                                stringsAsFactors = FALSE),
                     data.frame(region = regions[r], region_label = region_label[r],
                                process = "Colonization", year = years[-1],
                                mean = col_mean, lo = col_lo, hi = col_hi,
                                stringsAsFactors = FALSE),
                     data.frame(region = regions[r], region_label = region_label[r],
                                process = "Persistence", year = years[-1],
                                mean = per_mean, lo = per_lo, hi = per_hi,
                                stringsAsFactors = FALSE))
  }
  
  dyn_reg$process <- factor(dyn_reg$process, levels = c("Occupancy", "Colonization", "Persistence"))
  
  # --- regions where the species was never detected ---
  # The model fixes these at zero (assumed absent), so there is nothing
  # estimated to plot. Their lines are left blank and the panel is labelled
  # instead. The tables still contain the zeros.
  absent <- regions[n_det_reg == 0]
  
  occ_plot <- occ_reg
  occ_plot[occ_plot$region %in% absent, c("mean", "lo", "hi")] <- NA
  
  dyn_plot <- dyn_reg
  dyn_plot[dyn_plot$region %in% absent, c("mean", "lo", "hi")] <- NA
  
  # One text label per absent region (none if the species was detected everywhere)
  n_absent <- length(absent)
  absent_note <- data.frame(region_label = region_label[n_det_reg == 0],
                            process = factor(rep("Occupancy", n_absent), levels = levels(dyn_reg$process)),
                            year  = rep(mean(years), n_absent),
                            mean  = rep(0.5, n_absent),
                            label = rep("Not detected\n(assumed absent)", n_absent))
  
  # --- fig3: occupancy by region ---
  fig3 <- ggplot(occ_plot, aes(x = year, y = mean)) +
    geom_ribbon(data = occ_plot[!is.na(occ_plot$lo), ],
                aes(ymin = lo, ymax = hi), fill = LINE_COL, alpha = 0.2) +
    geom_line(colour = LINE_COL, linewidth = 1, na.rm = TRUE) +
    geom_point(colour = LINE_COL, size = 2, na.rm = TRUE) +
    geom_text(data = absent_note, aes(label = label), colour = "grey40", size = 3.5) +
    facet_wrap(~ region_label, ncol = 2) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = years) +
    labs(title = paste0(sp_title, ": occupancy by region"),
         subtitle = "Posterior mean (95% CI) across all sites in the region.",
         x = "Year", y = "Occupancy (\u03c8)") +
    theme_occ +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(sp_fig_dir, paste0("fig3_regional_occupancy_", SPECIES, ".png")),
         fig3, width = 9, height = 7, dpi = 300)
  
  # --- fig4: all three processes by region ---
  fig4 <- ggplot(dyn_plot, aes(x = year, y = mean, colour = process, fill = process)) +
    geom_ribbon(data = dyn_plot[!is.na(dyn_plot$lo), ],   # skip blank years
                aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_point(size = 1.5, na.rm = TRUE) +
    geom_text(data = absent_note, aes(label = label), colour = "grey40", size = 3.5) +
    facet_grid(process ~ region_label) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = years[seq(1, nyear, 2)]) +
    scale_colour_manual(values = c(Occupancy = "#1b9e77", Colonization = "#d95f02",
                                   Persistence = "#7570b3"), guide = "none") +
    scale_fill_manual(values = c(Occupancy = "#1b9e77", Colonization = "#d95f02",
                                 Persistence = "#7570b3"), guide = "none") +
    labs(title = paste0(sp_title, ": dynamics by region"),
         subtitle = paste0("Posterior mean (95% CI). Colonization and persistence are realised rates; persistence is blank where ",
                           "too few sites were occupied to estimate it."),
         x = "Year", y = "Probability") +
    theme_occ +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text.y = element_text(angle = 0))
  
  ggsave(file.path(sp_fig_dir, paste0("fig4_regional_dynamics_", SPECIES, ".png")),
         fig4, width = 12, height = 7.5, dpi = 300)
  
  # ============================================================================
  # 6. DETECTION BY YEAR
  # ============================================================================
  
  p_df <- data.frame(year = years,
                     mean = colMeans(sl$p.year),
                     lo = apply(sl$p.year, 2, quantile, 0.025),
                     hi = apply(sl$p.year, 2, quantile, 0.975))
  
  fig5 <- ggplot(p_df, aes(x = year, y = mean)) +
    geom_hline(yintercept = mean(sl$mean.p), linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, colour = "grey30") +
    geom_point(size = 2.5, colour = "grey10") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = years) +
    labs(title = paste0(sp_title, ": detection probability by year"),
         subtitle = "Per-visit detection at average covariate values (95% CI). Dashed line = overall mean.",
         x = "Year", y = "Detection probability (p)") +
    theme_occ
  
  ggsave(file.path(sp_fig_dir, paste0("fig5_detection_by_year_", SPECIES, ".png")),
         fig5, width = 7, height = 4.5, dpi = 300)
  
  # ============================================================================
  # 7. SITE-LEVEL OCCUPANCY
  # ============================================================================
  
  # --- p_occupied: posterior P(z = 1), conditional on the data ---
  p_occ <- colMeans(z_d)                       # [site, year]
  
  # --- psi_expected: model-predicted occupancy, NOT conditional on the site's
  #     own detections. Built by running the Markov recursion forward for every
  #     posterior draw:  psi[t] = psi[t-1] * phi[t-1] + (1 - psi[t-1]) * gamma[t-1]
  psi_exp <- array(NA_real_, dim = c(n_draws, nsite, nyear))
  
  lin_psi <- sl$a.psi[, d$region] +
    outer(sl$b.psi.elev,  d$elevation) +
    outer(sl$b.psi.water, d$dist_water)
  psi_exp[, , 1] <- plogis(lin_psi)
  
  for (t in 2:nyear) {
    lin_gam <- sl$a.gam[, d$region] + sl$eps.gam[, t - 1] +
      outer(sl$b.gam.elev,  d$elevation) +
      outer(sl$b.gam.water, d$dist_water) +
      outer(sl$b.gam.road,  d$dist_road) +
      outer(sl$b.gam.harv,  d$dist_harvest[, t - 1])
    lin_phi <- sl$a.phi[, d$region] + sl$eps.phi[, t - 1] +
      outer(sl$b.phi.elev,  d$elevation) +
      outer(sl$b.phi.water, d$dist_water) +
      outer(sl$b.phi.road,  d$dist_road) +
      outer(sl$b.phi.harv,  d$dist_harvest[, t - 1])
    psi_exp[, , t] <- psi_exp[, , t - 1] * plogis(lin_phi) +
      (1 - psi_exp[, , t - 1]) * plogis(lin_gam)
  }
  
  # Regions where the species was never detected are fixed at zero in the model
  psi_exp[, d$region %in% which(n_det_reg == 0), ] <- 0
  
  psi_exp_mean <- colMeans(psi_exp)
  psi_exp_lo   <- apply(psi_exp, c(2, 3), quantile, 0.025)
  psi_exp_hi   <- apply(psi_exp, c(2, 3), quantile, 0.975)
  
  # --- one row per site-year ---
  idx <- expand.grid(i = 1:nsite, t = 1:nyear)
  ij  <- cbind(idx$i, idx$t)
  
  site_year <- data.frame(
    species      = SPECIES,
    site         = d$sites[idx$i],
    GRTS_Cell_ID = site_info$GRTS_Cell_ID[idx$i],
    Quadrant     = site_info$Quadrant[idx$i],
    region       = site_info$region[idx$i],
    lat          = site_info$lat[idx$i],
    long         = site_info$long[idx$i],
    year         = years[idx$t],
    surveyed     = surveyed[ij],
    n_visits     = d$J[ij],
    detected     = ifelse(surveyed[ij], as.integer(detected[ij]), NA),
    p_occupied   = ifelse(surveyed[ij], round(p_occ[ij], 4), NA),
    psi_expected = round(psi_exp_mean[ij], 4),
    psi_expected_lo95 = round(psi_exp_lo[ij], 4),
    psi_expected_hi95 = round(psi_exp_hi[ij], 4),
    stringsAsFactors = FALSE
  )
  site_year <- site_year[order(site_year$site, site_year$year), ]
  
  # --- one row per site ---
  site_mean_draws <- rowMeans(psi_exp, dims = 2)   # [draws, site], averaged over years
  
  site_summary <- data.frame(
    species      = SPECIES,
    site         = d$sites,
    GRTS_Cell_ID = site_info$GRTS_Cell_ID,
    Quadrant     = site_info$Quadrant,
    region       = site_info$region,
    lat          = site_info$lat,
    long         = site_info$long,
    n_years_surveyed = rowSums(surveyed),
    n_years_detected = rowSums(detected),
    ever_detected    = ever_detected,
    mean_p_occupied  = round(rowSums(p_occ * surveyed) / rowSums(surveyed), 4),
    mean_psi_expected      = round(colMeans(site_mean_draws), 4),
    mean_psi_expected_lo95 = round(apply(site_mean_draws, 2, quantile, 0.025), 4),
    mean_psi_expected_hi95 = round(apply(site_mean_draws, 2, quantile, 0.975), 4),
    stringsAsFactors = FALSE
  )
  
  # --- fig6: map of mean P(occupied) across surveyed years ---
  map_df <- site_summary %>%
    mutate(status = ifelse(ever_detected, "Detected here", "Never detected")) %>%
    arrange(ever_detected)                     # detected sites drawn on top
  
  fig6 <- ggplot(map_df, aes(x = long, y = lat)) +
    geom_point(aes(fill = mean_p_occupied, shape = status),
               size = 3.2, colour = "grey25", stroke = 0.3, alpha = 0.9) +
    scale_shape_manual(values = c("Detected here" = 21, "Never detected" = 24), name = NULL) +
    scale_fill_distiller(palette = "OrRd", direction = 1, limits = c(0, 1),
                         name = "Mean\nP(occupied)") +
    coord_quickmap() +
    labs(title = paste0(sp_title, ": mean probability of occupancy per site, ",
                        min(years), "\u2013", max(years)),
         subtitle = "Averaged over surveyed years only. Colour on a triangle is model inference, not observation.",
         x = "Longitude", y = "Latitude") +
    theme_occ
  
  ggsave(file.path(sp_fig_dir, paste0("fig6_site_map_mean_", SPECIES, ".png")),
         fig6, width = 9, height = 8, dpi = 300)
  
  # --- fig7: map by year, surveyed site-years only ---
  map_year <- site_year %>%
    filter(surveyed) %>%
    mutate(status = ifelse(detected == 1, "Detected", "Not detected")) %>%
    arrange(year, detected)
  
  fig7 <- ggplot(map_year, aes(x = long, y = lat)) +
    geom_point(aes(fill = p_occupied, shape = status),
               size = 1.9, colour = "grey30", stroke = 0.2, alpha = 0.9) +
    facet_wrap(~ year, ncol = 3) +
    scale_shape_manual(values = c("Detected" = 21, "Not detected" = 24), name = NULL) +
    scale_fill_distiller(palette = "OrRd", direction = 1, limits = c(0, 1),
                         name = "P(occupied)") +
    coord_quickmap() +
    labs(title = paste0(sp_title, ": P(occupancy) by site and year"),
         subtitle = "Surveyed site-years only. Triangles = not detected that year; their colour is model inference.",
         x = NULL, y = NULL) +
    theme_occ +
    theme(legend.position = "bottom", axis.text = element_text(size = 7))
  
  ggsave(file.path(sp_fig_dir, paste0("fig7_site_map_by_year_", SPECIES, ".png")),
         fig7, width = 11, height = 12, dpi = 300)
  
  # ============================================================================
  # 8. TRENDS
  # ============================================================================
  # Slope = linear regression of occupancy on year, fitted to each posterior draw.
  # Change = last year minus first year. P(increase) = share of draws where
  # occupancy was higher in the last year than the first.
  
  year_c <- (1:nyear) - mean(1:nyear)
  
  trend_sources <- list(sl$psi.fs, sl$psi.surv, psi_regw)
  trend_names   <- c("BC (all sites)", "BC (surveyed sites)", "BC (region-weighted)")
  trend_ndet    <- c(n_det_total, n_det_total, n_det_total)
  for (r in 1:nregion) {
    trend_sources[[length(trend_sources) + 1]] <- sl$psi.reg[, r, ]
    trend_names <- c(trend_names, regions[r])
    trend_ndet  <- c(trend_ndet, n_det_reg[r])
  }
  
  trends <- data.frame()
  for (k in seq_along(trend_sources)) {
    m <- trend_sources[[k]]
    slope  <- as.vector(m %*% year_c) / sum(year_c^2)
    change <- m[, nyear] - m[, 1]
    trends <- rbind(trends, data.frame(
      species = SPECIES, scope = trend_names[k], n_detections = trend_ndet[k],
      psi_first = round(mean(m[, 1]), 4),
      psi_last  = round(mean(m[, nyear]), 4),
      slope = round(mean(slope), 5),
      slope_lo95 = round(quantile(slope, 0.025), 5),
      slope_hi95 = round(quantile(slope, 0.975), 5),
      change = round(mean(change), 4),
      change_lo95 = round(quantile(change, 0.025), 4),
      change_hi95 = round(quantile(change, 0.975), 4),
      p_increase = round(mean(change > 0), 4),
      stringsAsFactors = FALSE))
  }
  rownames(trends) <- NULL
  
  # ============================================================================
  # 9. EXPORTS
  # ============================================================================
  
  # Occupancy trajectories: three BC-wide versions plus each region
  traj <- rbind(
    data.frame(scope = "BC (all sites)", year = years,
               mean = colMeans(sl$psi.fs),
               lo95 = apply(sl$psi.fs, 2, quantile, 0.025),
               hi95 = apply(sl$psi.fs, 2, quantile, 0.975),
               naive = naive_bc, n_surveyed = colSums(surveyed)),
    data.frame(scope = "BC (surveyed sites)", year = years,
               mean = colMeans(sl$psi.surv),
               lo95 = apply(sl$psi.surv, 2, quantile, 0.025),
               hi95 = apply(sl$psi.surv, 2, quantile, 0.975),
               naive = naive_bc, n_surveyed = colSums(surveyed)),
    data.frame(scope = "BC (region-weighted)", year = years,
               mean = colMeans(psi_regw),
               lo95 = apply(psi_regw, 2, quantile, 0.025),
               hi95 = apply(psi_regw, 2, quantile, 0.975),
               naive = NA, n_surveyed = colSums(surveyed)),
    data.frame(scope = occ_reg$region, year = occ_reg$year,
               mean = occ_reg$mean, lo95 = occ_reg$lo, hi95 = occ_reg$hi,
               naive = occ_reg$naive, n_surveyed = occ_reg$n_surveyed)
  )
  traj$species <- SPECIES
  traj <- traj[, c("species", "scope", "year", "mean", "lo95", "hi95", "naive", "n_surveyed")]
  traj[, c("mean", "lo95", "hi95", "naive")] <- round(traj[, c("mean", "lo95", "hi95", "naive")], 4)
  
  # Colonisation and persistence
  dyn_out <- rbind(
    data.frame(scope = "BC (all sites)",
               process = as.character(dyn_bc$process), year = dyn_bc$year,
               mean = dyn_bc$mean, lo95 = dyn_bc$lo, hi95 = dyn_bc$hi),
    data.frame(scope = dyn_reg$region,
               process = as.character(dyn_reg$process), year = dyn_reg$year,
               mean = dyn_reg$mean, lo95 = dyn_reg$lo, hi95 = dyn_reg$hi)
  )
  dyn_out <- dyn_out[dyn_out$process != "Occupancy", ]
  dyn_out$species <- SPECIES
  dyn_out <- dyn_out[, c("species", "scope", "process", "year", "mean", "lo95", "hi95")]
  dyn_out[, c("mean", "lo95", "hi95")] <- round(dyn_out[, c("mean", "lo95", "hi95")], 4)
  
  cov_out <- cov_eff
  cov_out$species <- SPECIES
  cov_out$process <- as.character(cov_out$process)
  cov_out$covariate <- as.character(cov_out$covariate)
  cov_out <- cov_out[, c("species", "process", "covariate", "param",
                         "mean", "lo95", "hi95", "lo50", "hi50", "f", "excludes_zero")]
  cov_out[, c("mean", "lo95", "hi95", "lo50", "hi50", "f")] <-
    round(cov_out[, c("mean", "lo95", "hi95", "lo50", "hi50", "f")], 4)
  
  write.csv(cov_out,      file.path(per_sp_dir, paste0("covariate_effects_", SPECIES, ".csv")), row.names = FALSE)
  write.csv(traj,         file.path(per_sp_dir, paste0("occupancy_trajectory_", SPECIES, ".csv")), row.names = FALSE)
  write.csv(dyn_out,      file.path(per_sp_dir, paste0("dynamics_", SPECIES, ".csv")), row.names = FALSE)
  write.csv(trends,       file.path(per_sp_dir, paste0("trends_", SPECIES, ".csv")), row.names = FALSE)
  write.csv(site_year,    file.path(per_sp_dir, paste0("site_year_occupancy_", SPECIES, ".csv")), row.names = FALSE)
  write.csv(site_summary, file.path(per_sp_dir, paste0("site_summary_", SPECIES, ".csv")), row.names = FALSE)
  
  
  print(trends[, c("scope", "n_detections", "psi_first", "psi_last",
                   "slope", "slope_lo95", "slope_hi95", "p_increase")], row.names = FALSE)
  
  # free memory before the next species
  rm(fit, sl, z_d, psi_exp, psi_exp_lo, psi_exp_hi, site_mean_draws,
     lin_psi, lin_gam, lin_phi, col_bc, per_bc, col_r, per_r, psi_regw)
  gc()
}

# ============================================================================
# PART 2. COMBINE PER-SPECIES TABLES
# ============================================================================

export_types <- c("covariate_effects", "occupancy_trajectory", "dynamics",
                  "trends", "site_year_occupancy", "site_summary")

for (et in export_types) {
  files <- list.files(per_sp_dir, pattern = paste0("^", et, "_[A-Z]+\\.csv$"),
                      full.names = TRUE)
  if (length(files) == 0) next
  combined <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
  write.csv(combined, file.path(TAB_DIR, paste0("ALL_SPECIES_", et, ".csv")),
            row.names = FALSE)
  cat("Combined", length(files), "species ->", paste0("ALL_SPECIES_", et, ".csv"), "\n")
}

cov_all   <- read.csv(file.path(TAB_DIR, "ALL_SPECIES_covariate_effects.csv"), stringsAsFactors = FALSE)
trend_all <- read.csv(file.path(TAB_DIR, "ALL_SPECIES_trends.csv"), stringsAsFactors = FALSE)

# ============================================================================
# PART 3. CROSS-SPECIES FIGURES
# ============================================================================

# --- A. Covariate effects, all species ---
panel_order <- c("Initial occupancy: Elevation", "Initial occupancy: Distance to water",
                 "Colonization: Elevation", "Colonization: Distance to water",
                 "Colonization: Distance to road", "Colonization: Distance to harvest",
                 "Persistence: Elevation", "Persistence: Distance to water",
                 "Persistence: Distance to road", "Persistence: Distance to harvest",
                 "Detection: Clutter", "Detection: Temperature", "Detection: Julian date")

cov_all$panel <- factor(paste0(cov_all$process, ": ", cov_all$covariate), levels = panel_order)
cov_all$species <- factor(cov_all$species, levels = rev(sort(unique(cov_all$species))))
cov_all$direction <- factor(ifelse(cov_all$lo95 > 0, "Positive",
                                   ifelse(cov_all$hi95 < 0, "Negative", "Includes zero")),
                            levels = c("Positive", "Negative", "Includes zero"))

figA <- ggplot(cov_all, aes(x = mean, y = species, colour = direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_linerange(aes(xmin = lo95, xmax = hi95), linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_wrap(~ panel, ncol = 4, scales = "free_x") +
  scale_colour_manual(values = c(Positive = "#2166ac", Negative = "#b2182b",
                                 `Includes zero` = "grey70"),
                      labels = c(Positive = "Positive (95% CI above 0)",
                                 Negative = "Negative (95% CI below 0)",
                                 `Includes zero` = "95% CI includes 0"),
                      name = NULL, drop = FALSE) +
  labs(title = "Covariate effects across species (BC-wide)",
       subtitle = "Posterior mean and 95% credible interval. Logit scale, per 1 SD of covariate.",
       x = "Effect", y = NULL) +
  theme_occ +
  theme(legend.position = "bottom")

ggsave(file.path(all_fig_dir, "figA_covariate_effects_all_species.png"),
       figA, width = 13, height = 13, dpi = 300)

# --- B. BC-wide occupancy trend, all species ---
bc_trend <- trend_all %>%
  filter(scope == "BC (all sites)") %>%
  mutate(direction = case_when(p_increase >= 0.95 ~ "Increase (P \u2265 0.95)",
                               p_increase <= 0.05 ~ "Decrease (P \u2265 0.95)",
                               TRUE ~ "Uncertain"),
         species = reorder(species, slope))

figB <- ggplot(bc_trend, aes(x = slope, y = species, colour = direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_linerange(aes(xmin = slope_lo95, xmax = slope_hi95), linewidth = 0.6) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c("Increase (P \u2265 0.95)" = "#1b7837",
                                 "Decrease (P \u2265 0.95)" = "#b2182b",
                                 "Uncertain" = "grey50"), name = NULL) +
  labs(title = "BC-wide occupancy trend by species",
       subtitle = paste0("Linear slope of occupancy on year (95% CI), all sites. Colour: posterior probability that ",
                         "occupancy in the last year exceeds the first."),
       x = "Change in occupancy per year", y = NULL) +
  theme_occ +
  theme(legend.position = "bottom")

ggsave(file.path(all_fig_dir, "figB_bc_trends_all_species.png"),
       figB, width = 8, height = 6.5, dpi = 300)

# --- C. Regional occupancy trends, all species ---
# Species-region combinations with zero detections are dropped: there is no
# trend to estimate in a region where the species was never recorded.
reg_trend <- trend_all %>%
  filter(!grepl("^BC", scope), n_detections > 0) %>%
  mutate(direction = case_when(p_increase >= 0.95 ~ "Increase (P \u2265 0.95)",
                               p_increase <= 0.05 ~ "Decrease (P \u2265 0.95)",
                               TRUE ~ "Uncertain"))
reg_trend$species <- factor(reg_trend$species, levels = rev(sort(unique(reg_trend$species))))

figC <- ggplot(reg_trend, aes(x = slope, y = species, colour = direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_linerange(aes(xmin = slope_lo95, xmax = slope_hi95), linewidth = 0.6) +
  geom_point(size = 2) +
  facet_wrap(~ scope, nrow = 1) +
  scale_colour_manual(values = c("Increase (P \u2265 0.95)" = "#1b7837",
                                 "Decrease (P \u2265 0.95)" = "#b2182b",
                                 "Uncertain" = "grey50"), name = NULL) +
  labs(title = "Regional occupancy trends by species",
       subtitle = "Linear slope of occupancy on year (95% CI). Species with no detections in a region are omitted from that panel.",
       x = "Change in occupancy per year", y = NULL) +
  theme_occ +
  theme(legend.position = "bottom")

ggsave(file.path(all_fig_dir, "figC_regional_trends_all_species.png"),
       figC, width = 13, height = 6.5, dpi = 300)

# --- D. Occupancy in the final year, species x region ---
# A compact picture of where each species is and isn't.
last_year <- max(read.csv(file.path(TAB_DIR, "ALL_SPECIES_occupancy_trajectory.csv"))$year)
traj_all  <- read.csv(file.path(TAB_DIR, "ALL_SPECIES_occupancy_trajectory.csv"),
                      stringsAsFactors = FALSE)

heat <- traj_all %>%
  filter(year == last_year, !grepl("^BC", scope)) %>%
  left_join(trend_all %>% dplyr::select(species, scope, n_detections),
            by = c("species", "scope"))
heat$species <- factor(heat$species, levels = rev(sort(unique(heat$species))))

figD <- ggplot(heat, aes(x = scope, y = species, fill = mean)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = ifelse(n_detections == 0, "\u2013", sprintf("%.2f", mean))),
            size = 3, colour = ifelse(heat$mean > 0.55, "white", "grey15")) +
  scale_fill_distiller(palette = "OrRd", direction = 1, limits = c(0, 1),
                       name = paste0("Occupancy\n", last_year)) +
  labs(title = paste0("Estimated occupancy by region, ", last_year),
       subtitle = "Posterior mean across all sites in each region. \u2013 = species never detected in that region.",
       x = NULL, y = NULL) +
  theme_occ +
  theme(axis.line = element_blank(), axis.ticks = element_blank())

ggsave(file.path(all_fig_dir, "figD_occupancy_heatmap_all_species.png"),
       figD, width = 8, height = 6.5, dpi = 300)


