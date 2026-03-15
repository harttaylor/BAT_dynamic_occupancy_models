################################################################################
# 01a_prepare_data_BC.R — BC-Only Data Preparation
#
# Prepares species-level JAGS data for BC sites only.
# Includes distance-to-harvest as a covariate (available for all BC sites).
#
# Output: data/processed/BC/QUAD_[SPECIES]_data.rds
# Then run: 02_fit_models.R with INPUT_DIR = "data/processed/BC"
################################################################################

library(dplyr)
library(tidyr)

# ============================================================================
# SETTINGS
# ============================================================================

START_YEAR <- 2017
END_YEAR   <- 2025

RAW_DATA_FILE  <- "data/raw/final_detection_data_combined.csv"
DISTANCE_FILE  <- "data/BatSites_2025_distancetofootprint.csv"
ELEVATION_FILE <- "data/site_covariates_with_elevation.csv"
OUTPUT_DIR     <- "data/processed/BC"

VALID_QUADRANTS <- c("NE", "NW", "SE", "SW")
MAX_VISITS <- 4
MIN_VISITS <- 2
MIN_YEARS  <- 2

SPECIES_CODES <- c("ANPA", "COTO", "EPFU", "EUMA", "LABO", "LACI", "LANO",
                    "MYCA", "MYCI", "MYEV", "MYLU", "MYSE", "MYTH", "MYVO",
                    "MYYU", "PAHE", "TABR")

# ============================================================================
# 1. LOAD DATA
# ============================================================================

raw_data  <- read.csv(RAW_DATA_FILE)
site_covs <- read.csv(ELEVATION_FILE, stringsAsFactors = FALSE)

# Add distance covariates to BC if missing from elevation file
if (!"dist_road_km" %in% names(site_covs) || all(is.na(site_covs$dist_road_km))) {
  distance_covs <- read.csv(DISTANCE_FILE)
  distance_clean <- distance_covs %>%
    filter(Quadrant %in% VALID_QUADRANTS) %>%
    mutate(site = paste0(GRTS_Cell_, "_", Quadrant),
           dist_road_km = DIS_NearRd / 1000,
           dist_harvest_km = DIS_Harvst / 1000) %>%
    distinct(site, .keep_all = TRUE) %>%
    dplyr::select(site, dist_road_km, dist_harvest_km)
  site_covs <- left_join(site_covs, distance_clean, by = "site")
}

# ============================================================================
# 2. FILTER TO BC, VALID QUADRANTS, YEAR RANGE
# ============================================================================

bat_data <- raw_data %>%
  mutate(
    site   = paste0(GRTS.Cell.ID, "_", Quadrant),
    date   = as.Date(Night),
    year   = as.integer(format(date, "%Y")),
    julian = as.integer(format(date, "%j")),
    region = Region
  ) %>%
  filter(
    Quadrant %in% VALID_QUADRANTS,
    year >= START_YEAR, year <= END_YEAR,
    !is.na(region), region != "",
    region != "Alaska"
  )

# ============================================================================
# 3. EVERY-OTHER-NIGHT SUBSAMPLING
# ============================================================================

bat_data <- bat_data %>%
  group_by(site, year) %>% arrange(date) %>%
  mutate(night_num = row_number()) %>%
  filter(night_num %% 2 == 1) %>%
  mutate(visit = row_number()) %>%
  ungroup() %>%
  filter(visit <= MAX_VISITS) %>%
  dplyr::select(-night_num)

# ============================================================================
# 4. SITE FILTERING (min visits, min years)
# ============================================================================

site_year_visits <- bat_data %>%
  group_by(site, year) %>%
  summarize(n_visits = n(), .groups = "drop")
bat_data <- bat_data %>%
  semi_join(site_year_visits %>% filter(n_visits >= MIN_VISITS), by = c("site", "year"))

sites_keep <- bat_data %>%
  distinct(site, year) %>% count(site) %>%
  filter(n >= MIN_YEARS) %>% pull(site)
bat_data <- bat_data %>% filter(site %in% sites_keep)

all_sites <- sort(unique(bat_data$site))
all_years <- sort(unique(bat_data$year))
nsite <- length(all_sites); nyear <- length(all_years)

# ============================================================================
# 5. ALIGN & STANDARDIZE COVARIATES
# ============================================================================

site_covs <- site_covs %>%
  filter(site %in% all_sites) %>%
  arrange(match(site, all_sites))
stopifnot(nrow(site_covs) == nsite)
stopifnot(all(site_covs$site == all_sites))

# Impute NAs
site_covs <- site_covs %>% mutate(
  water_nearby    = ifelse(is.na(water_nearby), mean(water_nearby, na.rm = TRUE), water_nearby),
  dist_road_km    = ifelse(is.na(dist_road_km), mean(dist_road_km, na.rm = TRUE), dist_road_km),
  dist_harvest_km = ifelse(is.na(dist_harvest_km), mean(dist_harvest_km, na.rm = TRUE), dist_harvest_km),
  elevation_m     = ifelse(is.na(elevation_m), mean(elevation_m, na.rm = TRUE), elevation_m)
)

regions <- sort(unique(site_covs$region)); nregion <- length(regions)
site_covs$region_idx <- match(site_covs$region, regions)

# Standardize
safe_std <- function(x) {
  m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) s <- 1
  list(std = (x - m) / s, mean = m, sd = s)
}
elev_s <- safe_std(site_covs$elevation_m)
dist_road_s <- safe_std(site_covs$dist_road_km)
dist_harv_s <- safe_std(site_covs$dist_harvest_km)
water_mean <- mean(site_covs$water_nearby)
water_centered <- site_covs$water_nearby - water_mean

global_temp    <- mean(bat_data$Nightly.Mean.Temp, na.rm = TRUE)
global_clutter <- mean(bat_data$Percent.Clutter, na.rm = TRUE)
global_julian  <- mean(bat_data$julian, na.rm = TRUE)
if (is.na(global_temp)) global_temp <- 15
if (is.na(global_clutter)) global_clutter <- 50

std_params <- list(
  elevation_mean = elev_s$mean, elevation_sd = elev_s$sd,
  dist_road_mean = dist_road_s$mean, dist_road_sd = dist_road_s$sd,
  dist_harvest_mean = dist_harv_s$mean, dist_harvest_sd = dist_harv_s$sd,
  water_mean = water_mean,
  temp_mean = mean(bat_data$Nightly.Mean.Temp, na.rm = TRUE),
  temp_sd = sd(bat_data$Nightly.Mean.Temp, na.rm = TRUE),
  clutter_mean = mean(bat_data$Percent.Clutter, na.rm = TRUE),
  clutter_sd = sd(bat_data$Percent.Clutter, na.rm = TRUE),
  julian_mean = global_julian, julian_sd = sd(bat_data$julian, na.rm = TRUE)
)

# ============================================================================
# 6. BUILD & SAVE SPECIES DATA
# ============================================================================

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(std_params, file.path(OUTPUT_DIR, "standardization_params.rds"))

std_det <- function(arr, gmean, gsd) {
  if (is.na(gsd) || gsd == 0) gsd <- 1
  arr_std <- (arr - gmean) / gsd; arr_std[is.na(arr_std)] <- 0; arr_std
}

for (sp in SPECIES_CODES) {
  if (!sp %in% names(bat_data)) next
  y <- array(NA, dim = c(nsite, nyear, MAX_VISITS))
  temp_arr <- array(NA, dim = c(nsite, nyear, MAX_VISITS))
  clutter_arr <- array(NA, dim = c(nsite, nyear, MAX_VISITS))
  julian_arr <- array(NA, dim = c(nsite, nyear, MAX_VISITS))

  for (row in 1:nrow(bat_data)) {
    i <- match(bat_data$site[row], all_sites)
    k <- match(bat_data$year[row], all_years)
    j <- bat_data$visit[row]
    if (!is.na(i) && !is.na(k) && j <= MAX_VISITS) {
      sp_count <- bat_data[[sp]][row]
      y[i,k,j] <- ifelse(is.na(sp_count) || sp_count == 0, 0, 1)
      temp_arr[i,k,j] <- ifelse(is.na(bat_data$Nightly.Mean.Temp[row]),
                                 global_temp, bat_data$Nightly.Mean.Temp[row])
      clutter_arr[i,k,j] <- ifelse(is.na(bat_data$Percent.Clutter[row]),
                                    global_clutter, bat_data$Percent.Clutter[row])
      julian_arr[i,k,j] <- bat_data$julian[row]
    }
  }
  temp_arr    <- std_det(temp_arr, std_params$temp_mean, std_params$temp_sd)
  clutter_arr <- std_det(clutter_arr, std_params$clutter_mean, std_params$clutter_sd)
  julian_arr  <- std_det(julian_arr, std_params$julian_mean, std_params$julian_sd)

  J <- apply(y, 1:2, function(x) sum(!is.na(x)))
  nsurv <- array(NA, dim = c(nsite, nyear, MAX_VISITS))
  for (i in 1:nsite) for (k in 1:nyear) {
    valid <- which(!is.na(y[i,k,]))
    if (length(valid) > 0) nsurv[i,k,1:length(valid)] <- valid
  }

  n_det <- sum(y == 1, na.rm = TRUE)
  naive_occ <- mean(apply(y, 1, function(x) any(x == 1, na.rm = TRUE)), na.rm = TRUE)
  cat(sprintf("  %s: %d detections, naive occ = %.2f\n", sp, n_det, naive_occ))

  saveRDS(list(
    y = y, J = J, nsurv = nsurv,
    temp = temp_arr, clutter = clutter_arr, julian = julian_arr,
    nsite = nsite, nyear = nyear, nregion = nregion,
    elevation = elev_s$std, water = water_centered,
    dist_road = dist_road_s$std, dist_harvest = dist_harv_s$std,
    region = site_covs$region_idx,
    sites = all_sites, years = all_years, regions = regions,
    site_info = site_covs, std_params = std_params,
    psi_covariates = c("Intercept", "Elevation", "Water Nearby"),
    col_covariates = c("Intercept", "Dist to Road", "Dist to Harvest"),
    per_covariates = c("Intercept", "Dist to Road", "Dist to Harvest"),
    det_covariates = c("Intercept", "Clutter", "Temperature", "Julian")
  ), file.path(OUTPUT_DIR, paste0("QUAD_", sp, "_data.rds")))
}

write.csv(site_covs, file.path(OUTPUT_DIR, "site_covariates_final.csv"), row.names = FALSE)
cat("\n--- BC prep complete:", nsite, "sites x", nyear, "years ---\n")
cat("Output:", OUTPUT_DIR, "\n")
