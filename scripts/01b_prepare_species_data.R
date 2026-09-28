################################################################################
# 01b_prepare_species_data.R
#
# Builds the detection array for each species and packages it with the
# covariates from 01a into one .rds per species, ready for JAGS.
#
# Run 01a_prepare_covariates.R first.
#
################################################################################

library(dplyr)
library(tidyr)

# ============================================================================
# SETTINGS
# ============================================================================

# Two analyses use this script, one per report. Keep ONE block active and put
# a # in front of every line of the other. Use the same block in 01a, 01b, 02.

# --- BC report ---
DATA_DIR <- "data/processed/BC"

# --- Alaska report ---
# DATA_DIR <- "data/processed/BC_AK"

SPECIES_CODES <- c("ANPA", "COTO", "EPFU", "EUMA", "LABO", "LACI", "LANO",
                   "MYCA", "MYCI", "MYEV", "MYLU", "MYSE", "MYTH", "MYVO",
                   "MYYU", "PAHE", "TABR")

# Columns that are numeric and look like detections but are not species
NON_SPECIES_COLS <- c("F10K", "F20K", "F25K", "F30K", "F35K", "F40K",
                      "F45K", "F50K", "FeedingBuzz", "Flyingsquirrel",
                      "NoID", "NotBat", "Song", "Threebat", "Twobat", "nan")

# ============================================================================
# 1. LOAD THE ANALYSIS FRAME AND COVARIATES FROM 01a
# ============================================================================

bat_data  <- readRDS(file.path(DATA_DIR, "analysis_frame.rds"))
site_covs <- readRDS(file.path(DATA_DIR, "site_covariates.rds"))
covs      <- readRDS(file.path(DATA_DIR, "model_covariates.rds"))

all_sites  <- covs$all_sites
all_years  <- covs$all_years
nsite      <- covs$nsite
nyear      <- covs$nyear
nregion    <- covs$nregion
regions    <- covs$regions
MAX_VISITS <- covs$MAX_VISITS

# Detection covariate arrays: unsurveyed cells get 0 (the standardised mean).
# Those cells are never referenced by the likelihood; the fill just stops JAGS
# from seeing NAs in a data node.
temp_arr   <- covs$temp;   temp_arr[is.na(temp_arr)]     <- 0
julian_arr <- covs$julian; julian_arr[is.na(julian_arr)] <- 0

# ============================================================================
# 2. IDENTIFY AMBIGUOUS COLUMNS FOR EACH SPECIES
# ============================================================================
# An ambiguous column is any detection column whose name contains the species
# code but is not the bare code (LABOMYLU contains LABO; LABLPAHE does not).

all_cols <- names(bat_data)
candidate_cols <- all_cols[sapply(bat_data, is.numeric)]

not_species <- c(NON_SPECIES_COLS, "GRTS_Cell_ID", "Latitude", "Longitude",
                 "Elevation", "Distance_to_Clutter_m", "Percent_Clutter",
                 "Water_Nearby", "WaterDist", "Nightly_Mean_Temp",
                 "Nightly_Min_Temp", "Nightly_Max_Temp", "Nightly_Mean_RH",
                 "Nightly_Min_RH", "Nightly_Max_RH", "Nightly_Precipitation",
                 "Nightly_Mean_Windsp", "moon_illumination",
                 "moon_fraction_above_horizon", "dd", "DIST_WATER_M",
                 "DIST_ROAD", "DIST_HARVEST", "YEAR_HARVEST",
                 "year", "julian", "visit")

detection_cols <- setdiff(candidate_cols, not_species)

# ============================================================================
# 3. BUILD ONE DATASET PER SPECIES
# ============================================================================

summary_table <- data.frame()

for (sp in SPECIES_CODES) {
  
  if (!sp %in% names(bat_data)) next
  
  # ---- pure and ambiguous detection vectors, one value per row of the frame --
  pure <- bat_data[[sp]]
  pure[is.na(pure)] <- 0
  pure_det <- pure > 0
  
  amb_cols <- setdiff(detection_cols[grepl(sp, detection_cols, fixed = TRUE)], sp)
  
  if (length(amb_cols) > 0) {
    amb_mat <- as.matrix(bat_data[, amb_cols, drop = FALSE])
    amb_mat[is.na(amb_mat)] <- 0
    amb_det <- rowSums(amb_mat) > 0
  } else {
    amb_det <- rep(FALSE, nrow(bat_data))
  }
  
  amb_only <- amb_det & !pure_det
  
  # obs is 1 (detected) or 0 (not detected). Ambiguous-only visits stay 0.
  # Set obs[amb_only] <- NA_integer_ here to treat them as uninformative.
  obs <- as.integer(pure_det)
  
  # ---- fill the [site, year, visit] detection array ----
  y <- array(NA_integer_, dim = c(nsite, nyear, MAX_VISITS))
  
  i_idx <- match(bat_data$site, all_sites)
  k_idx <- match(bat_data$year, all_years)
  j_idx <- bat_data$visit
  
  y[cbind(i_idx, k_idx, j_idx)] <- obs
  
  # ---- J = number of informative visits per site-year ----
  # nsurv maps loop position j to the actual visit index.
  J <- apply(y, c(1, 2), function(x) sum(!is.na(x)))
  
  nsurv <- array(NA_integer_, dim = c(nsite, nyear, MAX_VISITS))
  for (i in 1:nsite) {
    for (k in 1:nyear) {
      valid <- which(!is.na(y[i, k, ]))
      if (length(valid) > 0) nsurv[i, k, 1:length(valid)] <- valid
    }
  }
  nsurv[is.na(nsurv)] <- 1   # placeholder, never referenced when J = 0
  
  # ---- surveyed indicator, for the surveyed-sites-only occupancy summary ----
  surv_ind <- (J > 0) * 1
  n_surv   <- colSums(surv_ind)
  
  # ---- region indicator matrix, for regional summaries ----
  reg_ind <- matrix(0, nrow = nsite, ncol = nregion)
  for (i in 1:nsite) reg_ind[i, covs$region_idx[i]] <- 1
  n_reg <- colSums(reg_ind)
  
  # ---- summary row ----
  det_reg   <- sapply(1:nregion, function(r) sum(y[covs$region_idx == r, , ] == 1, na.rm = TRUE))
  sites_reg <- sapply(1:nregion, function(r) {
    sum(apply(y[covs$region_idx == r, , , drop = FALSE], 1,
              function(x) any(x == 1, na.rm = TRUE)))
  })
  
  # Per-region counts become one column each. as.list() first, so the values
  # become columns rather than rows.
  region_tag <- gsub("[^A-Za-z0-9]+", "_", regions)
  det_cols   <- as.data.frame(as.list(setNames(det_reg,   paste0("det_",   region_tag))))
  site_cols  <- as.data.frame(as.list(setNames(sites_reg, paste0("sites_", region_tag))))
  
  base_row <- data.frame(
    species = sp,
    n_detections = sum(y == 1, na.rm = TRUE),
    n_informative_visits = sum(!is.na(y)),
    n_ambiguous_only = sum(amb_only),
    naive_occupancy = round(mean(apply(y, 1, function(x) any(x == 1, na.rm = TRUE))), 3),
    empty_regions = sum(det_reg == 0))
  
  summary_table <- rbind(summary_table, cbind(base_row, det_cols, site_cols))
  
  # ---- save ----
  saveRDS(list(
    species = sp,
    
    # detection data
    y = y, J = J, nsurv = nsurv,
    surv_ind = surv_ind, n_surv = n_surv,
    
    # dimensions
    nsite = nsite, nyear = nyear, nregion = nregion,
    sites = all_sites, years = all_years, regions = regions,
    
    # site covariates (standardised; dist_water is log(m + 1) then standardised)
    region = covs$region_idx,
    elevation = covs$elevation,
    dist_water = covs$dist_water,
    dist_road = covs$dist_road,
    clutter = covs$clutter,
    
    # site x year covariate
    dist_harvest = covs$dist_harvest,
    
    # detection covariates
    temp = temp_arr, julian = julian_arr,
    
    # regional bookkeeping
    reg_ind = reg_ind, n_reg = n_reg,
    
    # metadata
    site_info = site_covs,
    psi_covariates = c("Elevation", "Dist to Water"),
    col_covariates = c("Elevation", "Dist to Water", "Dist to Road", "Dist to Harvest"),
    per_covariates = c("Elevation", "Dist to Water", "Dist to Road", "Dist to Harvest"),
    det_covariates = c("Clutter", "Temperature", "Julian")
  ), file.path(DATA_DIR, paste0("QUAD_", sp, "_data.rds")))
}

write.csv(summary_table, file.path(DATA_DIR, "species_detection_summary.csv"),
          row.names = FALSE)

# ============================================================================
# 4. CHECK BEFORE MODELLING
# ============================================================================
# Species with <30 detections are reported descriptively, not modelled.
# Species with an empty region will fit, but that region's intercepts are
# prior-dominated and will show poor Rhat. Judge convergence there on the
# derived regional occupancy (psi.reg), not on a.psi/a.gam/a.phi.

summary_table[, c("species", "n_detections", "n_informative_visits",
                  "n_ambiguous_only", "naive_occupancy", "empty_regions")]
