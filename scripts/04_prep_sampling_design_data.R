################################################################################
# BC Bat Monitoring - Data Preparation: BC SITES ONLY (v5 - 2017 start)
# 
# READY TO RUN - Uses available covariates only
#
# Initial occupancy (ψ): Water Nearby (binary) + Distance to Road
# Colonization (γ): Distance to Road + Distance to Harvest  
# Persistence (φ): Distance to Road + Distance to Harvest
# Detection (p): Clutter + Temperature + Julian date
#
# Years: 2017-2024 (8 years, excludes unreliable 2016)
#
################################################################################

library(dplyr)
library(tidyr)

# ===== SETTINGS =====
OUTPUT_DIR <- "data/processed/bc_only_2017start"
OUTPUT_PREFIX <- "QUAD"
MIN_YEARS <- 2
MIN_VISITS <- 2
MAX_VISITS <- 4

START_YEAR <- 2017
END_YEAR <- 2024

VALID_QUADRANTS <- c("NE", "NW", "SE", "SW")

# ===== SITE CLASSIFICATION =====
ROTATING_GRTS <- c(
  # BC-WLRS
  186490, 137274, 77162, 243562, 226090, 128058, 114602, 143274,
  # BC Parks
  16426, 175354, 37274, 283586, 212906, 6202, 172090, 986, 
  265749, 317428, 38442, 132650, 115322, 206394, 26170, 
  230186, 17706, 145578, 58842, 110378
)

ALWAYS_ON_GRTS <- c(
  47210, 89962, 20330, 271731, 241002, 151210, 216746, 273050, 
  195242, 19514, 219626, 252394, 127850, 47322, 294832, 269777, 
  20730, 9466, 42026, 53210, 186330, 1194, 79018, 218282, 5994, 
  73002, 147306, 12522, 86074, 326072, 57322, 73786, 60586, 
  70250, 268521, 93290, 143338, 212266, 194666, 159722, 111466
)

cat("=" , rep("=", 60), "\n", sep = "")
cat("DATA PREPARATION: BC SITES ONLY (v5 - 2017 start)\n")
cat("Initial occupancy: Water (binary) + Distance to Road\n")
cat("Years:", START_YEAR, "-", END_YEAR, "\n")
cat("=" , rep("=", 60), "\n\n", sep = "")

# ===== 1. LOAD RAW DATA =====
raw_data <- read.csv("data/BatActivity_2016to2024.csv")

species_codes <- c("ANPA", "COTO", "EPFU", "EUMA", "LABO", "LACI", "LANO", 
                   "MYCA", "MYCI", "MYEV", "MYLU", "MYSE", "MYTH", "MYVO", 
                   "MYYU", "PAHE", "TABR")

# ===== 2. EXTRACT SITE-LEVEL COVARIATES FROM ALL RAW DATA =====
cat("Extracting site-level covariates from ALL raw data...\n")

raw_data <- raw_data %>%
  mutate(
    site = paste0(GRTS.Cell.ID, "_", Quadrant),
    quadrant = Quadrant,
    grts_cell = GRTS.Cell.ID,
    region = Region
  )

# Extract from ALL records (before filtering)
site_covariates_raw <- raw_data %>%
  filter(quadrant %in% VALID_QUADRANTS,
         !is.na(region),
         region != "Alaska",
         region != "") %>%
  group_by(site) %>%
  summarize(
    lat = mean(Lat, na.rm = TRUE),
    long = mean(Long, na.rm = TRUE),
    grts_cell = first(GRTS.Cell.ID),
    region = first(na.omit(region)),
    # Water nearby: binary (0/1) - take max across all records
    water_nearby = {
      w <- Water.Nearby[!is.na(Water.Nearby)]
      if(length(w) == 0) NA_real_ else as.numeric(max(w) > 0)
    },
    .groups = "drop"
  ) %>%
  mutate(
    site_type = case_when(
      grts_cell %in% ROTATING_GRTS ~ "Rotating",
      grts_cell %in% ALWAYS_ON_GRTS ~ "Always On",
      TRUE ~ "Always On"
    )
  )

cat("  Sites in raw data:", nrow(site_covariates_raw), "\n")
cat("  Sites with water data:", sum(!is.na(site_covariates_raw$water_nearby)), "\n")

# ===== 3. LOAD DISTANCE COVARIATES =====
cat("\nLoading distance covariates from footprint file...\n")

distance_covs <- read.csv("data/BatSites_2025_distancetofootprint.csv")

distance_clean <- distance_covs %>%
  filter(Quadrant %in% VALID_QUADRANTS) %>%
  mutate(
    site = paste0(GRTS_Cell_, "_", Quadrant),
    dist_road_km = DIS_NearRd / 1000,
    dist_harvest_km = DIS_Harvst / 1000
  ) %>%
  dplyr::select(site, dist_road_km, dist_harvest_km) %>%
  distinct(site, .keep_all = TRUE)

# Join to site info
site_covariates_raw <- site_covariates_raw %>%
  left_join(distance_clean, by = "site")

cat("  Sites with road distance:", sum(!is.na(site_covariates_raw$dist_road_km)), "\n")
cat("  Sites with harvest distance:", sum(!is.na(site_covariates_raw$dist_harvest_km)), "\n")

# ===== 4. FILTER DATA FOR ANALYSIS (2017-2024) =====
cat("\nFiltering data for analysis (2017-2024)...\n")


bat_data <- raw_data %>%
  mutate(
    date = as.Date(Night),
    year = as.integer(format(date, "%Y")),
    julian = as.integer(format(date, "%j"))
  ) %>%
  filter(year >= START_YEAR, year <= END_YEAR,
         quadrant %in% VALID_QUADRANTS,
         !is.na(region),
         region != "Alaska",
         region != "")

cat("  Records after filter:", nrow(bat_data), "\n")

# ===== 5. EVERY-OTHER-NIGHT SELECTION =====
bat_data <- bat_data %>%
  group_by(site, year) %>%
  arrange(date) %>%
  mutate(night_num = row_number()) %>%
  filter(night_num %% 2 == 1) %>%
  ungroup() %>%
  dplyr::select(-night_num)

# ===== 6. ASSIGN VISITS =====
bat_data <- bat_data %>%
  group_by(site, year) %>%
  arrange(date) %>%
  mutate(visit = row_number()) %>%
  ungroup() %>%
  filter(visit <= MAX_VISITS)

# ===== 7. FILTER SITES =====
site_year_visits <- bat_data %>%
  group_by(site, year) %>%
  summarize(n_visits = n(), .groups = "drop")

bat_data <- bat_data %>%
  semi_join(site_year_visits %>% filter(n_visits >= MIN_VISITS), 
            by = c("site", "year"))

site_years_count <- bat_data %>%
  distinct(site, year) %>%
  count(site)

sites_keep <- site_years_count %>%
  filter(n >= MIN_YEARS) %>%
  pull(site)

bat_data <- bat_data %>%
  filter(site %in% sites_keep)

# Master lists
all_sites <- sort(unique(bat_data$site))
all_years <- sort(unique(bat_data$year))
nsite <- length(all_sites)
nyear <- length(all_years)

cat("\nFinal: ", nsite, " sites x ", nyear, " years (", START_YEAR, "-", END_YEAR, ")\n", sep = "")

# ===== 8. SUBSET COVARIATES TO ANALYSIS SITES =====
site_covariates <- site_covariates_raw %>%
  filter(site %in% all_sites) %>%
  arrange(match(site, all_sites))

stopifnot(nrow(site_covariates) == nsite)
stopifnot(all(site_covariates$site == all_sites))

# ===== 9. REGION CODING =====
regions <- sort(unique(site_covariates$region))
nregion <- length(regions)

site_covariates <- site_covariates %>%
  mutate(region_idx = match(region, regions))

cat("Regions (", nregion, "): ", paste(regions, collapse = ", "), "\n", sep = "")

# ===== 10. COVARIATE COVERAGE =====
cat("\n*** COVARIATE COVERAGE ***\n")
cat("  Water nearby:", sum(!is.na(site_covariates$water_nearby)), "/", nsite, "\n")
cat("  Dist to road:", sum(!is.na(site_covariates$dist_road_km)), "/", nsite, "\n")
cat("  Dist to harvest:", sum(!is.na(site_covariates$dist_harvest_km)), "/", nsite, "\n")

# ===== 11. IMPUTE MISSING VALUES =====
water_global_mean <- mean(site_covariates$water_nearby, na.rm = TRUE)
dist_road_global_mean <- mean(site_covariates$dist_road_km, na.rm = TRUE)
dist_harvest_global_mean <- mean(site_covariates$dist_harvest_km, na.rm = TRUE)

if(is.na(water_global_mean)) water_global_mean <- 0.5
if(is.na(dist_road_global_mean)) dist_road_global_mean <- 1
if(is.na(dist_harvest_global_mean)) dist_harvest_global_mean <- 3

cat("\nImputing missing with global means:\n")
cat("  Water:", round(water_global_mean, 2), "\n")
cat("  Dist road:", round(dist_road_global_mean, 3), "km\n")
cat("  Dist harvest:", round(dist_harvest_global_mean, 2), "km\n")

site_covariates <- site_covariates %>%
  mutate(
    water_nearby = ifelse(is.na(water_nearby), water_global_mean, water_nearby),
    dist_road_km = ifelse(is.na(dist_road_km), dist_road_global_mean, dist_road_km),
    dist_harvest_km = ifelse(is.na(dist_harvest_km), dist_harvest_global_mean, dist_harvest_km)
  )

# ===== 12. STANDARDIZE COVARIATES =====
cat("\nStandardizing covariates...\n")

safe_standardize <- function(x, name = "variable") {
  x_mean <- mean(x, na.rm = TRUE)
  x_sd <- sd(x, na.rm = TRUE)
  if(is.na(x_sd) || x_sd == 0) {
    warning(paste(name, "has zero variance"))
    x_sd <- 1
  }
  x_std <- (x - x_mean) / x_sd
  x_std[!is.finite(x_std)] <- 0
  return(list(std = x_std, mean = x_mean, sd = x_sd))
}

# Water is binary - don't standardize, just center
water_mean <- mean(site_covariates$water_nearby)
water_centered <- site_covariates$water_nearby - water_mean

dist_road_scaled <- safe_standardize(site_covariates$dist_road_km, "dist_road")
dist_harvest_scaled <- safe_standardize(site_covariates$dist_harvest_km, "dist_harvest")

cat("  Dist road: mean=", round(dist_road_scaled$mean, 3), " SD=", round(dist_road_scaled$sd, 3), "\n", sep = "")
cat("  Dist harvest: mean=", round(dist_harvest_scaled$mean, 2), " SD=", round(dist_harvest_scaled$sd, 2), "\n", sep = "")

# Detection covariates
global_temp <- mean(bat_data$Nightly.Mean.Temp, na.rm = TRUE)
global_clutter <- mean(bat_data$Percent.Clutter, na.rm = TRUE)
if(is.na(global_temp)) global_temp <- 15
if(is.na(global_clutter)) global_clutter <- 50

# ===== 13. PREPARE SPECIES ARRAYS =====
cat("\nPreparing species arrays...\n")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

for(sp in species_codes) {
  
  y <- array(NA, dim = c(nsite, nyear, MAX_VISITS))
  n_det_covs <- 4
  x_p <- array(NA, dim = c(nsite, nyear, MAX_VISITS, n_det_covs))
  
  for(row in 1:nrow(bat_data)) {
    i <- match(bat_data$site[row], all_sites)
    k <- match(bat_data$year[row], all_years)
    j <- bat_data$visit[row]
    
    if(!is.na(i) && !is.na(k) && j <= MAX_VISITS) {
      sp_count <- bat_data[[sp]][row]
      y[i, k, j] <- ifelse(is.na(sp_count) || sp_count == 0, 0, 1)
      
      x_p[i, k, j, 1] <- 1
      x_p[i, k, j, 2] <- ifelse(is.na(bat_data$Percent.Clutter[row]), 
                                global_clutter, bat_data$Percent.Clutter[row])
      x_p[i, k, j, 3] <- ifelse(is.na(bat_data$Nightly.Mean.Temp[row]), 
                                global_temp, bat_data$Nightly.Mean.Temp[row])
      x_p[i, k, j, 4] <- bat_data$julian[row]
    }
  }
  
  J <- apply(y, 1:2, function(x) sum(!is.na(x)))
  
  nsurv <- array(NA, dim = c(nsite, nyear, MAX_VISITS))
  for(i in 1:nsite) {
    for(k in 1:nyear) {
      valid <- which(!is.na(y[i, k, ]))
      if(length(valid) > 0) nsurv[i, k, 1:length(valid)] <- valid
    }
  }
  
  # Standardize detection covariates
  for(cov in 2:n_det_covs) {
    vals <- x_p[,,,cov]
    vals_mean <- mean(vals, na.rm = TRUE)
    vals_sd <- sd(vals, na.rm = TRUE)
    if(is.na(vals_sd) || vals_sd == 0) vals_sd <- 1
    x_p[,,,cov] <- (vals - vals_mean) / vals_sd
  }
  x_p[is.na(x_p)] <- 0
  
  # Initial occupancy: intercept + water (centered) + dist_road (standardized)
  x_psi <- cbind(
    rep(1, nsite),
    water_centered,
    dist_road_scaled$std
  )
  
  # Colonization: intercept + dist_road + dist_harvest
  x_col <- cbind(
    rep(1, nsite),
    dist_road_scaled$std,
    dist_harvest_scaled$std
  )
  
  # Persistence: intercept + dist_road + dist_harvest
  x_per <- cbind(
    rep(1, nsite),
    dist_road_scaled$std,
    dist_harvest_scaled$std
  )
  
  n_det <- sum(y == 1, na.rm = TRUE)
  naive_occ <- mean(apply(y, 1:2, function(x) any(x == 1, na.rm = TRUE)), na.rm = TRUE)
  
  cat(sprintf("  %s: %d detections, naive occ = %.2f\n", sp, n_det, naive_occ))
  
  model_data <- list(
    y = y, J = J, nsurv = nsurv,
    nsite = nsite, nyear = nyear, nvisit = MAX_VISITS,
    x.psi = x_psi, nbeta.psi = ncol(x_psi),
    x.p = x_p, nbeta.p = n_det_covs,
    x.col = x_col, nbeta.col = ncol(x_col),
    x.per = x_per, nbeta.per = ncol(x_per),
    region = site_covariates$region_idx, 
    nregion = nregion,
    sites = all_sites, 
    years = all_years,
    site_info = site_covariates,
    regions = regions,
    scale = "QUADRANT",
    analysis_type = "BC_only_2017start_v5",
    start_year = START_YEAR,
    end_year = END_YEAR,
    rotating_grts = ROTATING_GRTS,
    always_on_grts = ALWAYS_ON_GRTS,
    psi_covariate_names = c("Intercept", "Water Nearby", "Dist to Road"),
    col_covariate_names = c("Intercept", "Dist to Road", "Dist to Harvest"),
    per_covariate_names = c("Intercept", "Dist to Road", "Dist to Harvest"),
    water_mean = water_mean,
    dist_road_mean = dist_road_scaled$mean, dist_road_sd = dist_road_scaled$sd,
    dist_harvest_mean = dist_harvest_scaled$mean, dist_harvest_sd = dist_harvest_scaled$sd
  )
  
  saveRDS(model_data, paste0(OUTPUT_DIR, "/", OUTPUT_PREFIX, "_", sp, "_model_data.rds"))
}

# Save site info
write.csv(site_covariates, paste0(OUTPUT_DIR, "/site_covariates_v5.csv"), row.names = FALSE)

cat("\n", rep("=", 60), "\n", sep = "")
cat("DATA PREPARATION COMPLETE (v5)\n")
cat("Sites:", nsite, "| Years:", nyear, "\n")
cat("Rotating:", sum(site_covariates$site_type == "Rotating"), 
    "| Always-on:", sum(site_covariates$site_type == "Always On"), "\n")
cat("\nCovariates:\n")
cat("  ψ (initial): Water Nearby (binary) + Dist to Road\n")
cat("  γ (colonization): Dist to Road + Dist to Harvest\n")
cat("  φ (persistence): Dist to Road + Dist to Harvest\n")
cat("\nOutput:", OUTPUT_DIR, "\n")
cat(rep("=", 60), "\n", sep = "")
