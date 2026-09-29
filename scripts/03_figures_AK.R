################################################################################
# 03_figures_AK.R
#
# Figures and tables for the Alaska sites, from the BC + Alaska
# run of 01a, 01b and 02. Only species detected in Alaska are plotted.
#
# PART 1 loops over species and makes the per-species figures and tables.
#
# PART 2 combines the per-species tables and makes the cross-species figures.
################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)

# ============================================================================
# SETTINGS
# ============================================================================

DATA_DIR <- "data/processed/BC_AK"
FIT_DIR  <- "outputs/BC_AK/fits"
OUT_DIR  <- "outputs/BC_AK"

FOCAL_REGION <- "Alaska"

det_summary     <- read.csv(file.path(DATA_DIR, "species_detection_summary.csv"))
SPECIES_TO_PLOT <- det_summary$species[det_summary$det_Alaska > 0]
# SPECIES_TO_PLOT <- "MYLU"      # uncomment to test one species

FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR <- file.path(OUT_DIR, "tables")
dir.create(file.path(FIG_DIR, "_all_species"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(TAB_DIR, "per_species"),  recursive = TRUE, showWarnings = FALSE)

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

cat("Species:", paste(SPECIES_TO_PLOT, collapse = ", "), "\n")

# ============================================================================
# PART 1. PER-SPECIES FIGURES AND TABLES
# ============================================================================

for (SPECIES in SPECIES_TO_PLOT) {
  
  # ---------------------------------------------------------------------------
  # 1.1 LOAD
  # ---------------------------------------------------------------------------
  
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
  sp_fig_dir <- file.path(FIG_DIR, SPECIES)
  dir.create(sp_fig_dir, recursive = TRUE, showWarnings = FALSE)
  
  nyear   <- d$nyear
  years   <- d$years
  n_draws <- length(sl$mean.p)
  
  # Alaska's region number and its sites. Everything below uses only these.
  r    <- match(FOCAL_REGION, d$regions)
  rows <- which(d$region == r)
  site_info <- d$site_info[rows, ]
  n_sites   <- length(rows)
  
  # ---------------------------------------------------------------------------
  # 1.2 OBSERVED DATA IN ALASKA
  # ---------------------------------------------------------------------------
  
  surveyed <- d$J[rows, , drop = FALSE] > 0                                    # [site, year]
  detected <- apply(d$y[rows, , , drop = FALSE], c(1, 2),
                    function(x) any(x == 1, na.rm = TRUE))                      # [site, year]
  detected[!surveyed] <- FALSE
  ever_detected <- apply(detected, 1, any)
  
  n_det      <- sum(d$y[rows, , ] == 1, na.rm = TRUE)
  n_surveyed <- colSums(surveyed)
  naive      <- ifelse(n_surveyed > 0, colSums(detected) / n_surveyed, NA)
  
  # A "detection night" is one survey night (visit) at one site with at least
  # one call identified to this species.
  region_note <- paste0(FOCAL_REGION, ": ", n_sites, " sites, detected at ", sum(ever_detected),
                        ", ", n_det, " detection nights")
  
  # ---------------------------------------------------------------------------
  # 1.3 FIG 1: COVARIATE EFFECTS
  # ---------------------------------------------------------------------------
  # The slopes are shared by every BC and Alaska site; only the intercepts differ
  # by region. "Alaska vs BC" under Detection is how much higher (positive) or
  # lower (negative) per-visit detection is at Alaska sites. Covariates are
  # standardised (distances log(x+1) first), so a slope is the change in
  # logit(probability) per 1 SD of the covariate.
  
  cov_eff <- data.frame(
    param = c("b.psi.elev", "b.psi.water",
              "b.gam.elev", "b.gam.water", "b.gam.road", "b.gam.harv",
              "b.phi.elev", "b.phi.water", "b.phi.road", "b.phi.harv",
              "b.p.ak", "b.p.clutter", "b.p.temp", "b.p.julian"),
    process = c(rep("Initial occupancy", 2), rep("Colonization", 4),
                rep("Persistence", 4), rep("Detection", 4)),
    covariate = c("Elevation", "Distance to water",
                  "Elevation", "Distance to water", "Distance to road", "Distance to harvest",
                  "Elevation", "Distance to water", "Distance to road", "Distance to harvest",
                  "Alaska vs BC", "Clutter", "Temperature", "Julian date"),
    stringsAsFactors = FALSE)
  
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
  cov_eff$direction <- factor(ifelse(cov_eff$lo95 > 0, "Positive",
                                     ifelse(cov_eff$hi95 < 0, "Negative", "Includes zero")),
                              levels = c("Positive", "Negative", "Includes zero"))
  cov_eff$process <- factor(cov_eff$process,
                            levels = c("Initial occupancy", "Colonization", "Persistence", "Detection"))
  cov_eff$covariate <- factor(cov_eff$covariate,
                              levels = rev(c("Elevation", "Distance to water", "Distance to road",
                                             "Distance to harvest", "Alaska vs BC", "Clutter",
                                             "Temperature", "Julian date")))
  
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
         subtitle = "Posterior mean and 95% credible interval. Slopes are shared by all BC and Alaska sites.",
         x = "Effect on logit scale (per 1 SD of covariate)", y = NULL) +
    theme_occ +
    theme(legend.position = "bottom", strip.text.y = element_text(angle = 0))
  
  ggsave(file.path(sp_fig_dir, paste0("fig1_covariate_effects_", SPECIES, ".png")),
         fig1, width = 8, height = 7, dpi = 300)
  
  # ---------------------------------------------------------------------------
  # 1.4 FIG 2: OCCUPANCY, COLONISATION AND PERSISTENCE IN ALASKA
  # ---------------------------------------------------------------------------
  # Colonisation and persistence are realised rates, calculated from the latent
  # occupancy states of the Alaska sites: the share of unoccupied sites that
  # became occupied, and of occupied sites that stayed occupied. Rates that are
  # undefined in most draws (e.g. persistence when almost no sites are
  # occupied) are left blank rather than estimated from noise.
  
  z_d   <- sl$z[, rows, , drop = FALSE]      # [draws, Alaska site, year]
  psi_r <- sl$psi.reg[, r, ]                 # [draws, year]
  
  col_r <- matrix(NA, n_draws, nyear - 1)
  per_r <- matrix(NA, n_draws, nyear - 1)
  for (t in 2:nyear) {
    prev <- z_d[, , t - 1, drop = FALSE]
    cur  <- z_d[, , t,     drop = FALSE]
    col_r[, t - 1] <- rowSums((1 - prev) * cur) / rowSums(1 - prev)
    per_r[, t - 1] <- rowSums(prev * cur) / rowSums(prev)
  }
  
  col_ok <- colMeans(is.finite(col_r)) >= 0.5
  per_ok <- colMeans(is.finite(per_r)) >= 0.5
  
  dyn <- rbind(
    data.frame(process = "Occupancy", year = years,
               mean = colMeans(psi_r),
               lo = apply(psi_r, 2, quantile, 0.025),
               hi = apply(psi_r, 2, quantile, 0.975),
               naive = naive),
    data.frame(process = "Colonization", year = years[-1],
               mean = ifelse(col_ok, colMeans(col_r, na.rm = TRUE), NA),
               lo = ifelse(col_ok, apply(col_r, 2, quantile, 0.025, na.rm = TRUE), NA),
               hi = ifelse(col_ok, apply(col_r, 2, quantile, 0.975, na.rm = TRUE), NA),
               naive = NA),
    data.frame(process = "Persistence", year = years[-1],
               mean = ifelse(per_ok, colMeans(per_r, na.rm = TRUE), NA),
               lo = ifelse(per_ok, apply(per_r, 2, quantile, 0.025, na.rm = TRUE), NA),
               hi = ifelse(per_ok, apply(per_r, 2, quantile, 0.975, na.rm = TRUE), NA),
               naive = NA))
  dyn$process <- factor(dyn$process, levels = c("Occupancy", "Colonization", "Persistence"))
  
  fig2 <- ggplot(dyn, aes(x = year, y = mean)) +
    geom_ribbon(data = dyn[!is.na(dyn$lo), ],             # skip blank years
                aes(ymin = lo, ymax = hi), fill = LINE_COL, alpha = 0.2) +
    geom_line(colour = LINE_COL, linewidth = 1, na.rm = TRUE) +
    geom_point(colour = LINE_COL, size = 2, na.rm = TRUE) +
    facet_wrap(~ process, nrow = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = years) +
    labs(title = paste0(sp_title, ": ", FOCAL_REGION),
         subtitle = paste0(region_note, ". Posterior mean (95% CI)."),
         x = "Year", y = "Probability") +
    theme_occ +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(sp_fig_dir, paste0("fig2_dynamics_", SPECIES, ".png")),
         fig2, width = 11, height = 4.5, dpi = 300)
  
  # ---------------------------------------------------------------------------
  # 1.5 FIG 3: DETECTION BY YEAR AT ALASKA SITES
  # ---------------------------------------------------------------------------
  # Alaska detection = BC detection + the Alaska offset (b.p.ak), on the logit
  # scale. The year-to-year pattern is shared with BC.
  
  p_ak <- plogis(sl$a.p + sl$eps.p + sl$b.p.ak)      # [draws, year]
  
  p_df <- data.frame(year = years,
                     mean = colMeans(p_ak),
                     lo = apply(p_ak, 2, quantile, 0.025),
                     hi = apply(p_ak, 2, quantile, 0.975))
  
  fig3 <- ggplot(p_df, aes(x = year, y = mean)) +
    geom_hline(yintercept = mean(plogis(sl$a.p + sl$b.p.ak)), linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, colour = "grey30") +
    geom_point(size = 2.5, colour = "grey10") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = years) +
    labs(title = paste0(sp_title, ": detection probability by year, ", FOCAL_REGION),
         subtitle = "Per-visit detection at Alaska sites, at average covariate values (95% CI). Dashed line = overall mean.",
         x = "Year", y = "Detection probability (p)") +
    theme_occ
  
  ggsave(file.path(sp_fig_dir, paste0("fig3_detection_by_year_", SPECIES, ".png")),
         fig3, width = 7, height = 4.5, dpi = 300)
  
  # ---------------------------------------------------------------------------
  # 1.6 SITE-LEVEL OCCUPANCY (Alaska sites only)
  # ---------------------------------------------------------------------------
  # p_occupied   = posterior P(site occupied that year), given its detections.
  # psi_expected = occupancy the model predicts from region, covariates and
  #                dynamics, NOT conditioned on the site's own detections.
  #                Built by running psi[t] = psi[t-1]*phi + (1 - psi[t-1])*gamma
  #                forward for every posterior draw.
  
  p_occ <- colMeans(z_d)                       # [site, year]
  
  reg   <- d$region[rows]
  elev  <- d$elevation[rows]
  water <- d$dist_water[rows]
  road  <- d$dist_road[rows]
  harv  <- d$dist_harvest[rows, , drop = FALSE]
  
  psi_exp <- array(NA_real_, dim = c(n_draws, n_sites, nyear))
  psi_exp[, , 1] <- plogis(sl$a.psi[, reg, drop = FALSE] +
                             outer(sl$b.psi.elev, elev) + outer(sl$b.psi.water, water))
  
  for (t in 2:nyear) {
    lin_gam <- sl$a.gam[, reg, drop = FALSE] + sl$eps.gam[, t - 1] +
      outer(sl$b.gam.elev, elev) + outer(sl$b.gam.water, water) +
      outer(sl$b.gam.road, road) + outer(sl$b.gam.harv, harv[, t - 1])
    lin_phi <- sl$a.phi[, reg, drop = FALSE] + sl$eps.phi[, t - 1] +
      outer(sl$b.phi.elev, elev) + outer(sl$b.phi.water, water) +
      outer(sl$b.phi.road, road) + outer(sl$b.phi.harv, harv[, t - 1])
    psi_exp[, , t] <- psi_exp[, , t - 1] * plogis(lin_phi) +
      (1 - psi_exp[, , t - 1]) * plogis(lin_gam)
  }
  
  psi_exp_mean <- colMeans(psi_exp)
  psi_exp_lo   <- apply(psi_exp, c(2, 3), quantile, 0.025)
  psi_exp_hi   <- apply(psi_exp, c(2, 3), quantile, 0.975)
  
  # --- one row per site-year ---
  idx <- expand.grid(i = 1:n_sites, t = 1:nyear)
  ij  <- cbind(idx$i, idx$t)
  
  site_year <- data.frame(
    species      = SPECIES,
    site         = d$sites[rows][idx$i],
    GRTS_Cell_ID = site_info$GRTS_Cell_ID[idx$i],
    Quadrant     = site_info$Quadrant[idx$i],
    region       = FOCAL_REGION,
    lat          = site_info$lat[idx$i],
    long         = site_info$long[idx$i],
    year         = years[idx$t],
    surveyed     = surveyed[ij],
    n_visits     = d$J[rows, , drop = FALSE][ij],
    detected     = ifelse(surveyed[ij], as.integer(detected[ij]), NA),
    p_occupied   = ifelse(surveyed[ij], round(p_occ[ij], 4), NA),
    psi_expected      = round(psi_exp_mean[ij], 4),
    psi_expected_lo95 = round(psi_exp_lo[ij], 4),
    psi_expected_hi95 = round(psi_exp_hi[ij], 4),
    stringsAsFactors = FALSE)
  site_year <- site_year[order(site_year$site, site_year$year), ]
  
  # --- one row per site ---
  site_mean_draws <- rowMeans(psi_exp, dims = 2)   # [draws, site], averaged over years
  
  site_summary <- data.frame(
    species      = SPECIES,
    site         = d$sites[rows],
    GRTS_Cell_ID = site_info$GRTS_Cell_ID,
    Quadrant     = site_info$Quadrant,
    region       = FOCAL_REGION,
    lat          = site_info$lat,
    long         = site_info$long,
    n_years_surveyed = rowSums(surveyed),
    n_years_detected = rowSums(detected),
    ever_detected    = ever_detected,
    mean_p_occupied  = round(rowSums(p_occ * surveyed) / rowSums(surveyed), 4),
    mean_psi_expected      = round(colMeans(site_mean_draws), 4),
    mean_psi_expected_lo95 = round(apply(site_mean_draws, 2, quantile, 0.025), 4),
    mean_psi_expected_hi95 = round(apply(site_mean_draws, 2, quantile, 0.975), 4),
    stringsAsFactors = FALSE)
  
  # ---------------------------------------------------------------------------
  # 1.7 FIG 4 AND 5: SITE MAPS
  # ---------------------------------------------------------------------------
  
  map_df <- site_summary %>%
    mutate(status = ifelse(ever_detected, "Detected here", "Never detected")) %>%
    arrange(ever_detected)                     # detected sites drawn on top
  
  fig4 <- ggplot(map_df, aes(x = long, y = lat)) +
    geom_point(aes(fill = mean_p_occupied, shape = status),
               size = 3.2, colour = "grey25", stroke = 0.3, alpha = 0.9) +
    scale_shape_manual(values = c("Detected here" = 21, "Never detected" = 24), name = NULL) +
    scale_fill_distiller(palette = "OrRd", direction = 1, limits = c(0, 1),
                         name = "Mean\nP(occupied)") +
    coord_quickmap() +
    labs(title = paste0(sp_title, ": mean probability of occupancy per site, ",
                        min(years), "\u2013", max(years)),
         subtitle = paste0(FOCAL_REGION, ". Averaged over surveyed years only. ",
                           "Colour on a triangle is model inference, not observation."),
         x = "Longitude", y = "Latitude") +
    theme_occ
  
  ggsave(file.path(sp_fig_dir, paste0("fig4_site_map_mean_", SPECIES, ".png")),
         fig4, width = 9, height = 8, dpi = 300)
  
  map_year <- site_year %>%
    filter(surveyed) %>%
    mutate(status = ifelse(detected == 1, "Detected", "Not detected")) %>%
    arrange(year, detected)
  
  fig5 <- ggplot(map_year, aes(x = long, y = lat)) +
    geom_point(aes(fill = p_occupied, shape = status),
               size = 2.2, colour = "grey30", stroke = 0.2, alpha = 0.9) +
    facet_wrap(~ year, ncol = 3) +
    scale_shape_manual(values = c("Detected" = 21, "Not detected" = 24), name = NULL) +
    scale_fill_distiller(palette = "OrRd", direction = 1, limits = c(0, 1),
                         name = "P(occupied)") +
    coord_quickmap() +
    labs(title = paste0(sp_title, ": P(occupancy) by site and year, ", FOCAL_REGION),
         subtitle = "Surveyed site-years only. Triangles = not detected that year; their colour is model inference.",
         x = NULL, y = NULL) +
    theme_occ +
    theme(legend.position = "bottom", axis.text = element_text(size = 7))
  
  ggsave(file.path(sp_fig_dir, paste0("fig5_site_map_by_year_", SPECIES, ".png")),
         fig5, width = 11, height = 12, dpi = 300)
  
  # ---------------------------------------------------------------------------
  # 1.8 TREND AND TABLES
  # ---------------------------------------------------------------------------
  # Slope = linear regression of occupancy on year, fitted to each posterior
  # draw. Change = last year minus first. P(increase) = share of draws where
  # occupancy was higher in the last year than the first.
  
  year_c <- (1:nyear) - mean(1:nyear)
  slope  <- as.vector(psi_r %*% year_c) / sum(year_c^2)
  change <- psi_r[, nyear] - psi_r[, 1]
  
  trend <- data.frame(
    species = SPECIES, region = FOCAL_REGION, n_sites = n_sites, n_detections = n_det,
    psi_first   = round(mean(psi_r[, 1]), 4),
    psi_last    = round(mean(psi_r[, nyear]), 4),
    slope       = round(mean(slope), 5),
    slope_lo95  = round(unname(quantile(slope, 0.025)), 5),
    slope_hi95  = round(unname(quantile(slope, 0.975)), 5),
    change      = round(mean(change), 4),
    change_lo95 = round(unname(quantile(change, 0.025)), 4),
    change_hi95 = round(unname(quantile(change, 0.975)), 4),
    p_increase  = round(mean(change > 0), 4))
  
  traj <- data.frame(species = SPECIES, region = FOCAL_REGION,
                     dyn[dyn$process == "Occupancy", c("year", "mean", "lo", "hi", "naive")],
                     n_surveyed = n_surveyed)
  names(traj)[names(traj) == "lo"] <- "lo95"
  names(traj)[names(traj) == "hi"] <- "hi95"
  
  dyn_out <- data.frame(species = SPECIES, region = FOCAL_REGION,
                        dyn[dyn$process != "Occupancy", c("process", "year", "mean", "lo", "hi")])
  names(dyn_out)[names(dyn_out) == "lo"] <- "lo95"
  names(dyn_out)[names(dyn_out) == "hi"] <- "hi95"
  
  cov_out <- data.frame(species = SPECIES,
                        process = as.character(cov_eff$process),
                        covariate = as.character(cov_eff$covariate),
                        cov_eff[, c("param", "mean", "lo95", "hi95", "lo50", "hi50", "f", "excludes_zero")])
  
  traj[, c("mean", "lo95", "hi95", "naive")] <- round(traj[, c("mean", "lo95", "hi95", "naive")], 4)
  dyn_out[, c("mean", "lo95", "hi95")]       <- round(dyn_out[, c("mean", "lo95", "hi95")], 4)
  cov_out[, c("mean", "lo95", "hi95", "lo50", "hi50", "f")] <-
    round(cov_out[, c("mean", "lo95", "hi95", "lo50", "hi50", "f")], 4)
  
  per_sp <- file.path(TAB_DIR, "per_species")
  write.csv(cov_out,      file.path(per_sp, paste0("covariate_effects_",    SPECIES, ".csv")), row.names = FALSE)
  write.csv(traj,         file.path(per_sp, paste0("occupancy_trajectory_", SPECIES, ".csv")), row.names = FALSE)
  write.csv(dyn_out,      file.path(per_sp, paste0("dynamics_",             SPECIES, ".csv")), row.names = FALSE)
  write.csv(trend,        file.path(per_sp, paste0("trends_",               SPECIES, ".csv")), row.names = FALSE)
  write.csv(site_year,    file.path(per_sp, paste0("site_year_occupancy_",  SPECIES, ".csv")), row.names = FALSE)
  write.csv(site_summary, file.path(per_sp, paste0("site_summary_",         SPECIES, ".csv")), row.names = FALSE)
  
  print(trend[, c("region", "n_detections", "psi_first", "psi_last",
                  "slope", "slope_lo95", "slope_hi95", "p_increase")], row.names = FALSE)
  
  # free memory before the next species
  rm(fit, sl, z_d, psi_exp, psi_exp_lo, psi_exp_hi, site_mean_draws, lin_gam, lin_phi)
  gc()
}

# ============================================================================
# PART 2. COMBINE TABLES AND MAKE CROSS-SPECIES FIGURES
# ============================================================================

per_sp <- file.path(TAB_DIR, "per_species")

for (tab in c("covariate_effects", "occupancy_trajectory", "dynamics",
              "trends", "site_year_occupancy", "site_summary")) {
  files <- list.files(per_sp, pattern = paste0("^", tab, "_[A-Z]+\\.csv$"), full.names = TRUE)
  combined <- do.call(rbind, lapply(files, read.csv))
  write.csv(combined, file.path(TAB_DIR, paste0("ALL_SPECIES_", tab, ".csv")), row.names = FALSE)
}

cov_all   <- read.csv(file.path(TAB_DIR, "ALL_SPECIES_covariate_effects.csv"))
trend_all <- read.csv(file.path(TAB_DIR, "ALL_SPECIES_trends.csv"))
traj_all  <- read.csv(file.path(TAB_DIR, "ALL_SPECIES_occupancy_trajectory.csv"))

fig_all <- file.path(FIG_DIR, "_all_species")

# --- A. Covariate effects, all species ---
panel_order <- c("Initial occupancy: Elevation", "Initial occupancy: Distance to water",
                 "Colonization: Elevation", "Colonization: Distance to water",
                 "Colonization: Distance to road", "Colonization: Distance to harvest",
                 "Persistence: Elevation", "Persistence: Distance to water",
                 "Persistence: Distance to road", "Persistence: Distance to harvest",
                 "Detection: Alaska vs BC", "Detection: Clutter", "Detection: Temperature",
                 "Detection: Julian date")

cov_all$panel     <- factor(paste0(cov_all$process, ": ", cov_all$covariate), levels = panel_order)
cov_all$species   <- factor(cov_all$species, levels = rev(sort(unique(cov_all$species))))
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
  labs(title = "Covariate effects, species detected in Alaska",
       subtitle = "Posterior mean and 95% credible interval. Logit scale, per 1 SD. Slopes shared by all BC and Alaska sites.",
       x = "Effect", y = NULL) +
  theme_occ +
  theme(legend.position = "bottom", strip.text = element_text(size = 9))

ggsave(file.path(fig_all, "figA_covariate_effects_all_species.png"),
       figA, width = 13, height = 11, dpi = 300)

# --- B. Occupancy trend in Alaska, all species ---
trend_all$direction <- ifelse(trend_all$p_increase >= 0.95, "Increase (P \u2265 0.95)",
                              ifelse(trend_all$p_increase <= 0.05, "Decrease (P \u2265 0.95)",
                                     "Uncertain"))
trend_all$species <- reorder(trend_all$species, trend_all$slope)

figB <- ggplot(trend_all, aes(x = slope, y = species, colour = direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_linerange(aes(xmin = slope_lo95, xmax = slope_hi95), linewidth = 0.6) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c("Increase (P \u2265 0.95)" = "#1b7837",
                                 "Decrease (P \u2265 0.95)" = "#b2182b",
                                 "Uncertain" = "grey50"), name = NULL) +
  labs(title = "Occupancy trend by species, Alaska",
       subtitle = "Linear slope of occupancy on year (95% CI). Colour: posterior probability that occupancy in the last year exceeds the first.",
       x = "Change in occupancy per year", y = NULL) +
  theme_occ +
  theme(legend.position = "bottom")

ggsave(file.path(fig_all, "figB_trends_all_species.png"),
       figB, width = 8, height = 5, dpi = 300)

# --- C. Occupancy trajectory in Alaska, all species ---
figC <- ggplot(traj_all, aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = LINE_COL, alpha = 0.2) +
  geom_line(colour = LINE_COL, linewidth = 1) +
  facet_wrap(~ species) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Occupancy by species, Alaska",
       subtitle = "Posterior mean (95% CI).",
       x = "Year", y = "Occupancy (\u03c8)") +
  theme_occ +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_all, "figC_occupancy_all_species.png"),
       figC, width = 10, height = 7, dpi = 300)


