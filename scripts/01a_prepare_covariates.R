################################################################################
# 01a_prepare_covariates.R
#
# COVARIATE PREPARATION ONLY. This script defines the analysis frame (which
# site-years and which visits are used) and prepares every covariate the models
# need. Doesn't touch species detections - that gets prepared in 01b
#
# So when covariates change (new years of harvest data, new landscape variable, 
# corrected clutter values) only this script needs editing, and 01b/02 don't need 
# modification to still work/run
#
# For the Alaska report, Alaska is simply a fifth region. The only Alaska-
# specific lines are the AK_ site ID prefix (section 2) and filling missing
# temperatures within each region (section 9)
#
# WHAT THIS PRODUCES (in data/processed/BC/ or data/processed/BC_AK/):
#   analysis_frame.rds        night-level rows that survive filtering, with
#                             site / year / visit already assigned
#   site_covariates.rds       one row per site, standardised + raw
#   siteyear_covariates.rds   site x year covariates (harvest), gap-filled
#   detection_covariates.rds  temp and julian arrays [site, year, visit]
#   standardization_params.rds  means and SDs, so new data can be put on the
#                             same scale later
#   site_covariates_final.csv human-readable version for collaborators
#
# COVARIATE DECISIONS MADE HERE (and why):
#   - DIST_ROAD, DIST_HARVEST, DIST_WATER_M are log(x+1) transformed before
#     standardising. All three are severely right-skewed (harvest reaches
#     79 km against a median of ~1.2 km). On the raw scale, standardisation is
#     driven by a handful of extreme sites and the fitted slope describes those
#     outliers rather than the bulk of the data.
#   - Water is DIST_WATER_M (continuous, 0% missing), replacing the old binary
#     Water_Nearby, which was missing at 38 sites in every year.
#   - WATER_TYPE and ROAD_TYPE are not used.
#   - Percent_Clutter is missing at 40 sites in EVERY year, so it cannot be
#     carried forward. Those sites get the mean of all sites and are flagged in
#     clutter_imputed so the decision stays visible.
#   - Nightly_Mean_Temp is 21% missing and genuinely varies by night, so it is
#     filled from the site-year mean, then the region-year mean, then the
#     region mean (Alaska is much cooler than BC, so a BC-wide fill would bias it).
#   - DIST_HARVEST is the only site covariate that varies by year, so it is stored as a site x year matrix. Years with no survey get
#     the nearest known value for that site carried forward, then backward.
#   - dd (degree days) is night-level and collinear with Julian date, so it is
#     not used.
#   - Habitat_Type is free text with 31 inconsistent levels, so it is not used.
################################################################################

library(dplyr)
library(tidyr)

# ============================================================================
# SETTINGS
# ============================================================================

# Two analyses use this script, one per report. Keep ONE block active and put
# a # in front of every line of the other. Use the same block in 01a, 01b, 02.

# --- BC report: BC sites only ---
RAW_DATA_FILE   <- "data/raw/AllYearActivitybyNight_2016_2025_BC_withWeather_NEAR.csv"
OUTPUT_DIR      <- "data/processed/BC"
EXCLUDE_REGIONS <- c("Alaska")

# --- Alaska report: BC and Alaska together (run 00_harmonise_BC_AK.R first) ---
# RAW_DATA_FILE   <- "data/processed/AllYearActivitybyNight_BC_AK_combined.csv"
# OUTPUT_DIR      <- "data/processed/BC_AK"
# EXCLUDE_REGIONS <- c()

START_YEAR <- 2017    # 2016 excluded: sparse and unreliable
END_YEAR   <- 2025

VALID_QUADRANTS <- c("NE", "NW", "SE", "SW")   # exact codes only (NE1, SW2 etc. are dropped)

MAX_VISITS <- 4       # every-other-night gives 4 visits from a 7-8 night deployment
MIN_VISITS <- 2       # a site-year needs 2 visits to inform detection
MIN_YEARS  <- 2       # a site needs 2 surveyed years to inform the dynamics

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# 1. LOAD AND NORMALISE COLUMN NAMES
# ============================================================================
# Header conventions have varied between exports (spaces, dots, underscores).
# Normalising here means nothing downstream has to guess.

raw_data <- read.csv(RAW_DATA_FILE, check.names = FALSE, stringsAsFactors = FALSE)

names(raw_data) <- gsub("[^A-Za-z0-9]+", "_", names(raw_data))
names(raw_data) <- gsub("_+$", "", names(raw_data))

cat("Loaded", nrow(raw_data), "rows,", ncol(raw_data), "columns\n")

# ============================================================================
# 2. BUILD SITE / YEAR / DATE FIELDS
# ============================================================================

bat_data <- raw_data

bat_data$site   <- paste0(bat_data$GRTS_Cell_ID, "_", bat_data$Quadrant)
# BC and Alaska number their GRTS cells separately, so the same number can
# appear in both. Alaska sites get an AK_ prefix (e.g. AK_206_NE) so they are
# never merged with a BC site. BC site IDs are unchanged.
is_ak <- bat_data$Region %in% "Alaska"
bat_data$site[is_ak] <- paste0("AK_", bat_data$site[is_ak])
bat_data$date   <- as.Date(bat_data$Night)
bat_data$year   <- as.integer(format(bat_data$date, "%Y"))
bat_data$julian <- as.integer(format(bat_data$date, "%j"))
bat_data$region <- bat_data$Region

n_bad_date <- sum(is.na(bat_data$date))
n_odd_year <- sum(bat_data$year < 2016 | bat_data$year > END_YEAR, na.rm = TRUE)
cat("Unparseable dates:", n_bad_date, "| rows outside 2016-", END_YEAR, ":", n_odd_year, "\n")

# ============================================================================
# 3. FILTER TO THE ANALYSIS FRAME
# ============================================================================

bat_data <- bat_data %>%
  filter(
    Quadrant %in% VALID_QUADRANTS,
    !is.na(year), year >= START_YEAR, year <= END_YEAR,
    !is.na(region), region != "",
    !region %in% EXCLUDE_REGIONS
  )

cat("After region / quadrant / year filter:", nrow(bat_data), "rows\n")

# ============================================================================
# 4. EVERY-OTHER-NIGHT SUBSAMPLING AND VISIT CAP
# ============================================================================
# Deployments run 7-8 consecutive nights. Taking every other night keeps the
# repeat visits approximately independent (NABat convention). Capping at 4 also
# means the visits always fall within roughly the first week of a deployment,
# which keeps the within-year closure assumption tight even at the 109
# site-years whose full deployment spans more than 14 days.

bat_data <- bat_data %>%
  group_by(site, year) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(night_num = row_number()) %>%
  filter(night_num %% 2 == 1) %>%
  mutate(visit = row_number()) %>%
  ungroup() %>%
  filter(visit <= MAX_VISITS) %>%
  dplyr::select(-night_num)

cat("After subsampling and visit cap:", nrow(bat_data), "rows\n")

# ============================================================================
# 5. DROP THIN SITE-YEARS AND THIN SITES
# ============================================================================

site_year_visits <- bat_data %>%
  group_by(site, year) %>%
  summarize(n_visits = n(), .groups = "drop")

bat_data <- bat_data %>%
  semi_join(site_year_visits %>% filter(n_visits >= MIN_VISITS),
            by = c("site", "year"))

sites_keep <- bat_data %>%
  distinct(site, year) %>%
  count(site) %>%
  filter(n >= MIN_YEARS) %>%
  pull(site)

bat_data <- bat_data %>% filter(site %in% sites_keep)

# Visit numbers must be contiguous after the filters above
bat_data <- bat_data %>%
  group_by(site, year) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(visit = row_number()) %>%
  ungroup()

all_sites <- sort(unique(bat_data$site))
all_years <- sort(unique(bat_data$year))
nsite <- length(all_sites)
nyear <- length(all_years)

cat("\nFINAL ANALYSIS FRAME\n")
cat("  Sites:", nsite, "\n")
cat("  Years:", nyear, "(", min(all_years), "-", max(all_years), ")\n")
cat("  Site-years surveyed:", nrow(distinct(bat_data, site, year)), "of", nsite * nyear, "\n")
cat("  Visits:", nrow(bat_data), "\n")

# ============================================================================
# 6. SITE-LEVEL COVARIATES
# ============================================================================
# These do not vary across years (confirmed in 00_data_diagnostics.R), so take
# the first non-missing value per site.

first_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  x[1]
}

site_covs <- bat_data %>%
  group_by(site) %>%
  summarize(
    GRTS_Cell_ID    = first(GRTS_Cell_ID),
    Quadrant        = first(Quadrant),
    region          = first(region),
    lat             = first_value(Latitude),
    long            = first_value(Longitude),
    elevation_m     = first_value(Elevation),
    dist_water_m    = first_value(DIST_WATER_M),
    dist_road_m     = first_value(DIST_ROAD),
    pct_clutter_raw = first_value(Percent_Clutter),
    .groups = "drop"
  ) %>%
  arrange(match(site, all_sites))

stopifnot(nrow(site_covs) == nsite)
stopifnot(all(site_covs$site == all_sites))

# --- clutter: missing at every year for some sites, so mean-impute and flag ---
site_covs$clutter_imputed <- is.na(site_covs$pct_clutter_raw)
clutter_mean_raw <- mean(site_covs$pct_clutter_raw, na.rm = TRUE)
site_covs$pct_clutter <- ifelse(site_covs$clutter_imputed,
                                clutter_mean_raw,
                                site_covs$pct_clutter_raw)

cat("\nClutter imputed at", sum(site_covs$clutter_imputed), "of", nsite,
    "sites (mean of all sites =", round(clutter_mean_raw, 1), ")\n")

# --- confirm nothing else is missing ---
cat("\nMissing values in site covariates:\n")
print(colSums(is.na(site_covs[, c("lat", "long", "elevation_m", "dist_water_m",
                                  "dist_road_m", "pct_clutter")])))

# --- log-transform the skewed distance covariates ---
site_covs$log_dist_water <- log(site_covs$dist_water_m + 1)
site_covs$log_dist_road  <- log(site_covs$dist_road_m  + 1)

# --- region index ---
regions <- sort(unique(site_covs$region))
nregion <- length(regions)
site_covs$region_idx <- match(site_covs$region, regions)

cat("\nSites per region:\n")
print(table(site_covs$region))

# ============================================================================
# 7. SITE x YEAR COVARIATES (harvest distance)
# ============================================================================
# DIST_HARVEST changes across years at some sites as new cutblocks appear.
# The model needs a value for EVERY site-year, including unsurveyed ones,
# because colonisation and persistence are defined for all years. Gaps are
# filled by carrying the nearest known value forward, then backward.

harvest_obs <- bat_data %>%
  group_by(site, year) %>%
  summarize(dist_harvest_m = first_value(DIST_HARVEST),
            year_harvest   = first_value(YEAR_HARVEST),
            .groups = "drop")

harvest_grid <- expand.grid(site = all_sites, year = all_years,
                            stringsAsFactors = FALSE) %>%
  left_join(harvest_obs, by = c("site", "year")) %>%
  arrange(site, year) %>%
  group_by(site) %>%
  fill(dist_harvest_m, year_harvest, .direction = "downup") %>%
  ungroup()

cat("\nHarvest site-years filled by carry-forward/backward:",
    sum(is.na(harvest_obs$dist_harvest_m[match(
      paste(harvest_grid$site, harvest_grid$year),
      paste(harvest_obs$site, harvest_obs$year))])), "\n")
cat("Harvest values still missing after fill:", sum(is.na(harvest_grid$dist_harvest_m)), "\n")

harvest_grid$log_dist_harvest <- log(harvest_grid$dist_harvest_m + 1)

# ============================================================================
# 8. STANDARDISE
# ============================================================================
# Site covariates are standardised across sites. Harvest is standardised across
# all site-years pooled, so a single slope is comparable across years.

safe_std <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) s <- 1
  list(std = (x - m) / s, mean = m, sd = s)
}

elev_s    <- safe_std(site_covs$elevation_m)
water_s   <- safe_std(site_covs$log_dist_water)
road_s    <- safe_std(site_covs$log_dist_road)
clutter_s <- safe_std(site_covs$pct_clutter)
harv_s    <- safe_std(harvest_grid$log_dist_harvest)

site_covs$elevation_std    <- elev_s$std
site_covs$dist_water_std   <- water_s$std
site_covs$dist_road_std    <- road_s$std
site_covs$clutter_std      <- clutter_s$std
harvest_grid$dist_harvest_std <- harv_s$std

# Reshape harvest to a [site, year] matrix in analysis-frame order
dist_harvest_mat <- harvest_grid %>%
  dplyr::select(site, year, dist_harvest_std) %>%
  pivot_wider(names_from = year, values_from = dist_harvest_std) %>%
  arrange(match(site, all_sites))

stopifnot(all(dist_harvest_mat$site == all_sites))
dist_harvest_mat <- as.matrix(dist_harvest_mat[, as.character(all_years)])
rownames(dist_harvest_mat) <- all_sites

stopifnot(!any(is.na(dist_harvest_mat)))

# ============================================================================
# 9. NIGHT-LEVEL DETECTION COVARIATES
# ============================================================================
# Temperature is missing on some nights. Fill from the site-year mean (nights
# within a site-year are within about a week of each other), then the mean for
# that region and year, then the region mean, then the overall mean.
# Julian date has no missing values.

temp_filled <- bat_data %>%
  group_by(site, year) %>%
  mutate(temp_f = ifelse(is.na(Nightly_Mean_Temp),
                         mean(Nightly_Mean_Temp, na.rm = TRUE),
                         Nightly_Mean_Temp)) %>%
  ungroup() %>%
  group_by(region, year) %>%
  mutate(temp_f = ifelse(is.na(temp_f),
                         mean(Nightly_Mean_Temp, na.rm = TRUE),
                         temp_f)) %>%
  ungroup() %>%
  group_by(region) %>%
  mutate(temp_f = ifelse(is.na(temp_f),
                         mean(Nightly_Mean_Temp, na.rm = TRUE),
                         temp_f)) %>%
  ungroup()

global_temp <- mean(bat_data$Nightly_Mean_Temp, na.rm = TRUE)
temp_filled$temp_f[is.na(temp_filled$temp_f) | is.nan(temp_filled$temp_f)] <- global_temp

cat("\nTemperature: ", sum(is.na(bat_data$Nightly_Mean_Temp)), " of ", nrow(bat_data),
    " nights filled (", round(100 * mean(is.na(bat_data$Nightly_Mean_Temp)), 1),
    "%)\n", sep = "")

temp_s   <- safe_std(temp_filled$temp_f)
julian_s <- safe_std(temp_filled$julian)

temp_filled$temp_std   <- temp_s$std
temp_filled$julian_std <- julian_s$std

# --- build [site, year, visit] arrays ---
temp_arr   <- array(NA_real_, dim = c(nsite, nyear, MAX_VISITS))
julian_arr <- array(NA_real_, dim = c(nsite, nyear, MAX_VISITS))

i_idx <- match(temp_filled$site, all_sites)
k_idx <- match(temp_filled$year, all_years)
j_idx <- temp_filled$visit

temp_arr[cbind(i_idx, k_idx, j_idx)]   <- temp_filled$temp_std
julian_arr[cbind(i_idx, k_idx, j_idx)] <- temp_filled$julian_std

# Cells with no survey stay NA here; 01b replaces them with 0 (the standardised
# mean) before they reach JAGS. They are never referenced by the likelihood.

# ============================================================================
# 10. SAVE
# ============================================================================

std_params <- list(
  elevation_mean = elev_s$mean,    elevation_sd = elev_s$sd,
  log_water_mean = water_s$mean,   log_water_sd = water_s$sd,
  log_road_mean  = road_s$mean,    log_road_sd  = road_s$sd,
  log_harvest_mean = harv_s$mean,  log_harvest_sd = harv_s$sd,
  clutter_mean   = clutter_s$mean, clutter_sd   = clutter_s$sd,
  temp_mean      = temp_s$mean,    temp_sd      = temp_s$sd,
  julian_mean    = julian_s$mean,  julian_sd    = julian_s$sd,
  clutter_mean_raw = clutter_mean_raw,
  note = "Distances are log(x+1) transformed BEFORE standardising."
)

saveRDS(bat_data,         file.path(OUTPUT_DIR, "analysis_frame.rds"))
saveRDS(site_covs,        file.path(OUTPUT_DIR, "site_covariates.rds"))
saveRDS(harvest_grid,     file.path(OUTPUT_DIR, "siteyear_covariates.rds"))
saveRDS(std_params,       file.path(OUTPUT_DIR, "standardization_params.rds"))

saveRDS(list(
  all_sites = all_sites, all_years = all_years,
  nsite = nsite, nyear = nyear, nregion = nregion, regions = regions,
  MAX_VISITS = MAX_VISITS,
  region_idx      = site_covs$region_idx,
  elevation       = site_covs$elevation_std,
  dist_water      = site_covs$dist_water_std,
  dist_road       = site_covs$dist_road_std,
  clutter         = site_covs$clutter_std,
  dist_harvest    = dist_harvest_mat,
  temp            = temp_arr,
  julian          = julian_arr
), file.path(OUTPUT_DIR, "model_covariates.rds"))

write.csv(site_covs, file.path(OUTPUT_DIR, "site_covariates_final.csv"),
          row.names = FALSE)

cat("\n--- Covariate prep complete ---\n")
cat("Output:", OUTPUT_DIR, "\n")
cat("Next: 01b_prepare_species_data.R\n")

# ============================================================================
# 11. QUICK SANITY CHECKS
# ============================================================================

cat("\nStandardised covariate summaries (should be mean ~0, sd ~1):\n")
print(round(data.frame(
  elevation    = c(mean(site_covs$elevation_std),  sd(site_covs$elevation_std)),
  dist_water   = c(mean(site_covs$dist_water_std), sd(site_covs$dist_water_std)),
  dist_road    = c(mean(site_covs$dist_road_std),  sd(site_covs$dist_road_std)),
  clutter      = c(mean(site_covs$clutter_std),    sd(site_covs$clutter_std)),
  dist_harvest = c(mean(dist_harvest_mat),         sd(dist_harvest_mat)),
  row.names = c("mean", "sd")
), 3))

cat("\nCorrelation among site covariates (watch for collinearity with region):\n")
print(round(cor(site_covs[, c("elevation_std", "dist_water_std",
                              "dist_road_std", "clutter_std")]), 2))

cat("\nMean standardised elevation by region (large spread = region/elevation confounding):\n")
print(round(tapply(site_covs$elevation_std, site_covs$region, mean), 2))