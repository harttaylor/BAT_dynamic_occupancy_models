################################################################################
# BC Bat Monitoring - Data Preparation for Dynamic Occupancy (JAGS)
# 
# SITE DEFINITION: GRTS.Cell.ID_Quadrant (core quadrants only: NE, NW, SE, SW)
# 
# COVARIATE STRUCTURE:
#   - Initial occupancy (psi): water [i]
#   - Col/Pers (gamma/phi): clutter + water [i, k] (site-year level)
#   - Detection (p): clutter + temp + julian + moon [i, k, j] (survey level)
#
################################################################################

library(dplyr)
library(tidyr)
library(purrr)

# ===== 1. LOAD AND CLEAN DATA =====

raw_data <- read.csv("data/BatActivity_2016to2024.csv")

species_codes <- c("ANPA", "COTO", "EPFU", "EUMA", "LABO", "LACI", "LANO", 
                   "MYCA", "MYCI", "MYEV", "MYLU", "MYSE", "MYTH", "MYVO", 
                   "MYYU", "PAHE", "TABR")

# Parse dates and create site identifier
raw_data <- raw_data %>%
  mutate(
    date = as.Date(Night),
    year = as.integer(format(date, "%Y")),
    julian = as.integer(format(date, "%j")),
    site = paste0(GRTS.Cell.ID, "_", Quadrant),
    park = SiteName,
    region = Region
  ) %>%
  filter(year >= 2016, year <= 2024) %>%
  arrange(site, year, date)

cat("=== DATA OVERVIEW ===\n")
cat("Years: 2016-2024\n")
cat("Total rows:", nrow(raw_data), "\n")
cat("Unique sites (all quadrants):", n_distinct(raw_data$site), "\n\n")

# ===== 2. FILTER TO CORE QUADRANTS (NE, NW, SE, SW) =====

cat("=== FILTERING TO CORE QUADRANTS ===\n")

core_quadrants <- c("NE", "NW", "SE", "SW")
bat_data <- raw_data %>%
  filter(Quadrant %in% core_quadrants)

cat("Keeping quadrants:", paste(core_quadrants, collapse = ", "), "\n")
cat("Rows after filter:", nrow(bat_data), "(dropped", nrow(raw_data) - nrow(bat_data), ")\n")
cat("Unique sites:", n_distinct(bat_data$site), "\n\n")

# ===== 3. ASSIGN VISITS WITHIN EACH SITE-YEAR =====
# Each night at a site is a visit (row = site x night)

cat("=== ASSIGNING VISITS ===\n")

# ------------------------------------------------------------------------------
# OPTION A: USE 4 visits per year
# ------------------------------------------------------------------------------

MAX_VISITS <- 4  # Most site-years have ~7 nights

bat_data <- bat_data %>%
  group_by(site, year) %>%
  arrange(date) %>%
  mutate(visit = row_number()) %>%
  ungroup()

# Cap at MAX_VISITS
bat_data <- bat_data %>%
  filter(visit <= MAX_VISITS)

# ------------------------------------------------------------------------------
# OPTION B: USE EVERY 2ND NIGHT (REDUCES AUTOCORRELATION)
# ------------------------------------------------------------------------------

# MAX_VISITS <- 4  # With every 2nd night from ~7 nights, expect ~3-4 visits
# 
# bat_data <- bat_data %>%
#   group_by(site, year) %>%
#   arrange(date) %>%
#   mutate(night_num = row_number()) %>%
#   # Keep every 2nd night (1st, 3rd, 5th, 7th)
#   filter(night_num %% 2 == 1) %>%
#   mutate(visit = row_number()) %>%
#   ungroup() %>%
#   select(-night_num)
# 
# # Cap at MAX_VISITS
# bat_data <- bat_data %>%
#   filter(visit <= MAX_VISITS)

# ------------------------------------------------------------------------------

# Check visit distribution
visit_dist <- bat_data %>%
  group_by(site, year) %>%
  summarize(n_visits = max(visit), .groups = "drop") %>%
  count(n_visits)

cat("\nVisits per site-year distribution:\n")
print(visit_dist)

# ===== 4. FILTER SITES WITH SUFFICIENT DATA =====

cat("\n=== FILTERING SITES ===\n")

# Remove site-years with < 2 visits
min_visits_required <- 2
site_year_visits <- bat_data %>%
  group_by(site, year) %>%
  summarize(n_visits = n(), .groups = "drop")

site_years_keep <- site_year_visits %>%
  filter(n_visits >= min_visits_required)

bat_data <- bat_data %>%
  semi_join(site_years_keep, by = c("site", "year"))

cat("Retained site-years with >=", min_visits_required, "visits:", nrow(site_years_keep), "\n")

# Remove sites with < 2 years of data (needed for colonization/extinction)
site_years_count <- bat_data %>%
  distinct(site, year) %>%
  count(site)

sites_keep <- site_years_count %>%
  filter(n >= 2) %>%
  pull(site)

bat_data <- bat_data %>%
  filter(site %in% sites_keep)

cat("Retained sites with >= 2 years:", length(sites_keep), "\n")

# Create master lists
all_sites <- sort(unique(bat_data$site))
all_years <- sort(unique(bat_data$year))
nsite <- length(all_sites)
nyear <- length(all_years)

cat("\nFinal dimensions:", nsite, "sites x", nyear, "years x", MAX_VISITS, "max visits\n")

# Final visit distribution
final_visit_dist <- bat_data %>%
  group_by(site, year) %>%
  summarize(n_visits = n(), .groups = "drop") %>%
  count(n_visits)

cat("\nFinal visits per site-year:\n")
print(final_visit_dist)

# ===== 5. SITE-LEVEL COVARIATES (for initial occupancy) =====
# Use ANY available value from any year for the site

cat("\n=== SITE-LEVEL COVARIATES ===\n")

site_covariates <- bat_data %>%
  group_by(site) %>%
  summarize(
    park = first(na.omit(park)),
    region = first(na.omit(region)),
    lat = mean(Lat, na.rm = TRUE),
    long = mean(Long, na.rm = TRUE),
    # For water: take the most common non-NA value
    water_raw = {
      water_vals <- Water.Nearby[!is.na(Water.Nearby)]
      if(length(water_vals) == 0) NA_real_ else {
        as.numeric(names(sort(table(water_vals), decreasing = TRUE)[1]))
      }
    },
    .groups = "drop"
  ) %>%
  mutate(
    water = ifelse(is.na(water_raw) | water_raw < 0, 0, water_raw)
  ) %>%
  arrange(match(site, all_sites))

cat("Sites with water = 1:", sum(site_covariates$water > 0), "\n")
cat("Sites with water = 0 (incl. imputed):", sum(site_covariates$water == 0), "\n")

# ===== 6. SITE-YEAR COVARIATES (for colonization/persistence) =====
# Hierarchy: site-year value -> site-level value (any year) -> global mean

cat("\n=== SITE-YEAR COVARIATES ===\n")

# Get site-year level means
site_year_raw <- bat_data %>%
  group_by(site, year) %>%
  summarize(
    clutter_mean = mean(Percent.Clutter, na.rm = TRUE),
    water_mean = {
      water_vals <- Water.Nearby[!is.na(Water.Nearby)]
      if(length(water_vals) == 0) NA_real_ else max(water_vals)
    },
    .groups = "drop"
  ) %>%
  mutate(
    clutter_mean = ifelse(is.nan(clutter_mean), NA, clutter_mean)
  )

# Get site-level values (from any year) for imputation
site_level_covs <- bat_data %>%
  group_by(site) %>%
  summarize(
    clutter_site = mean(Percent.Clutter, na.rm = TRUE),
    water_site = {
      water_vals <- Water.Nearby[!is.na(Water.Nearby)]
      if(length(water_vals) == 0) NA_real_ else max(water_vals)
    },
    .groups = "drop"
  ) %>%
  mutate(
    clutter_site = ifelse(is.nan(clutter_site), NA, clutter_site)
  )

# Global means for final fallback
global_clutter <- mean(bat_data$Percent.Clutter, na.rm = TRUE)
global_water <- 0

# Create complete site-year grid
site_year_grid <- expand.grid(
  site = all_sites,
  year = all_years,
  stringsAsFactors = FALSE
) %>%
  as_tibble()

# Merge and fill
site_year_covs <- site_year_grid %>%
  left_join(site_year_raw, by = c("site", "year")) %>%
  left_join(site_level_covs, by = "site") %>%
  mutate(
    # Use site-year value if available, else site-level, else global
    clutter_filled = case_when(
      !is.na(clutter_mean) ~ clutter_mean,
      !is.na(clutter_site) ~ clutter_site,
      TRUE ~ global_clutter
    ),
    water_filled = case_when(
      !is.na(water_mean) ~ water_mean,
      !is.na(water_site) ~ water_site,
      TRUE ~ global_water
    )
  )

n_sy_clutter <- sum(!is.na(site_year_covs$clutter_mean))
n_site_clutter <- sum(is.na(site_year_covs$clutter_mean) & !is.na(site_year_covs$clutter_site))
n_global_clutter <- sum(is.na(site_year_covs$clutter_mean) & is.na(site_year_covs$clutter_site))

cat("Clutter imputation:\n")
cat("  Site-year value:", n_sy_clutter, "\n")
cat("  Site-level (from other years):", n_site_clutter, "\n")
cat("  Global mean:", n_global_clutter, "\n")

# ===== 7. PREPARE ARRAYS FOR EACH SPECIES =====

cat("\n=== PREPARING JAGS DATA ===\n\n")

# Global means for detection covariate imputation
global_temp <- mean(bat_data$Nightly.Mean.Temp, na.rm = TRUE)
global_julian <- mean(bat_data$julian, na.rm = TRUE)
global_moon <- mean(bat_data$mphase, na.rm = TRUE)
global_clutter_det <- mean(bat_data$Percent.Clutter, na.rm = TRUE)

for(sp in species_codes) {
  
  # ----- Detection array y[i, k, j] -----
  y <- array(NA, dim = c(nsite, nyear, MAX_VISITS))
  
  # Detection covariates x.p[i, k, j, covariate]
  n_det_covs <- 5  # intercept, clutter, temp, julian, moon
  x_p <- array(NA, dim = c(nsite, nyear, MAX_VISITS, n_det_covs))
  
  det_impute_count <- 0
  
  # Fill in detection data and covariates
  for(row in 1:nrow(bat_data)) {
    site_name <- bat_data$site[row]
    year_val <- bat_data$year[row]
    visit_num <- bat_data$visit[row]
    
    i <- match(site_name, all_sites)
    k <- match(year_val, all_years)
    j <- visit_num
    
    if(!is.na(i) && !is.na(k) && j <= MAX_VISITS) {
      # Detection: 1 if species count > 0, else 0
      species_count <- bat_data[[sp]][row]
      y[i, k, j] <- ifelse(is.na(species_count), 0, ifelse(species_count > 0, 1, 0))
      
      # Detection covariates
      x_p[i, k, j, 1] <- 1  # Intercept
      
      # Clutter
      clutter_val <- bat_data$Percent.Clutter[row]
      if(!is.na(clutter_val)) {
        x_p[i, k, j, 2] <- clutter_val
      } else {
        x_p[i, k, j, 2] <- global_clutter_det
        det_impute_count <- det_impute_count + 1
      }
      
      # Temperature
      temp_val <- bat_data$Nightly.Mean.Temp[row]
      if(!is.na(temp_val)) {
        x_p[i, k, j, 3] <- temp_val
      } else {
        x_p[i, k, j, 3] <- global_temp
        det_impute_count <- det_impute_count + 1
      }
      
      # Julian day
      x_p[i, k, j, 4] <- bat_data$julian[row]
      
      # Moon phase
      moon_val <- bat_data$mphase[row]
      if(!is.na(moon_val)) {
        x_p[i, k, j, 5] <- moon_val
      } else {
        x_p[i, k, j, 5] <- global_moon
        det_impute_count <- det_impute_count + 1
      }
    }
  }
  
  # ----- Survey indices -----
  J <- apply(y, 1:2, function(x) sum(!is.na(x)))
  
  nsurv <- array(NA, dim = c(nsite, nyear, MAX_VISITS))
  for(i in 1:nsite) {
    for(k in 1:nyear) {
      valid <- which(!is.na(y[i, k, ]))
      if(length(valid) > 0) nsurv[i, k, 1:length(valid)] <- valid
    }
  }
  
  # ----- Indicator matrix (any detection in site-year) -----
  ind <- apply(y, c(1, 2), function(x) {
    if(all(is.na(x))) return(0)
    return(ifelse(any(x == 1, na.rm = TRUE), 1, 0))
  })
  
  # ----- Initial occupancy covariates x.psi[i, covariate] -----
  x_psi <- cbind(
    rep(1, nsite),
    site_covariates$water
  )
  
  # ----- Col/Pers covariates x.phi[i, k, covariate] and x.gamma[i, k, covariate] -----
  n_colpers_covs <- 3  # intercept, clutter, water
  x_phi <- array(NA, dim = c(nsite, nyear, n_colpers_covs))
  x_gamma <- array(NA, dim = c(nsite, nyear, n_colpers_covs))
  
  for(i in 1:nsite) {
    site_name <- all_sites[i]
    for(k in 1:nyear) {
      year_val <- all_years[k]
      
      sy_row <- site_year_covs %>% 
        filter(site == site_name, year == year_val)
      
      x_phi[i, k, 1] <- 1
      x_phi[i, k, 2] <- sy_row$clutter_filled
      x_phi[i, k, 3] <- sy_row$water_filled
      
      x_gamma[i, k, 1] <- 1
      x_gamma[i, k, 2] <- sy_row$clutter_filled
      x_gamma[i, k, 3] <- sy_row$water_filled
    }
  }
  
  # ----- Standardize continuous covariates -----
  
  # Detection covariates (columns 2-5)
  for(cov_idx in 2:n_det_covs) {
    vals <- x_p[, , , cov_idx]
    sd_val <- sd(vals, na.rm = TRUE)
    if(!is.na(sd_val) && sd_val > 0) {
      x_p[, , , cov_idx] <- (vals - mean(vals, na.rm = TRUE)) / sd_val
    }
  }
  
  # Col/pers covariates - standardize clutter (column 2)
  clutter_vals <- x_phi[, , 2]
  sd_clutter <- sd(clutter_vals, na.rm = TRUE)
  if(!is.na(sd_clutter) && sd_clutter > 0) {
    clutter_std <- (clutter_vals - mean(clutter_vals, na.rm = TRUE)) / sd_clutter
    x_phi[, , 2] <- clutter_std
    x_gamma[, , 2] <- clutter_std
  }
  
  # ----- Diagnostics -----
  n_surveyed <- sum(!is.na(y))
  n_detections <- sum(y == 1, na.rm = TRUE)
  n_zeros <- sum(y == 0, na.rm = TRUE)
  n_na <- sum(is.na(y))
  pct_surveyed <- round(100 * n_surveyed / length(y), 1)
  
  cat(sprintf("%-6s: %5d det, %6d non-det, %5d NA (%.1f%% surveyed)\n",
              sp, n_detections, n_zeros, n_na, pct_surveyed))
  
  # ----- Save model data -----
  
  model_data <- list(
    # Detection data
    y = y,
    nsurv = nsurv,
    J = J,
    ind = ind,
    
    # Dimensions
    nsite = nsite,
    nyear = nyear,
    nvisit = MAX_VISITS,
    
    # Covariate arrays
    x.psi = x_psi,
    nbeta.psi = ncol(x_psi),
    
    x.phi = x_phi,
    nbeta.phi = dim(x_phi)[3],
    
    x.gamma = x_gamma,
    nbeta.gamma = dim(x_gamma)[3],
    
    x.p = x_p,
    nbeta.p = dim(x_p)[4],
    
    # Metadata
    sites = all_sites,
    years = all_years,
    site_info = site_covariates,
    det_imputed = det_impute_count
  )
  
  # Create output directory if needed
  dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
  saveRDS(model_data, paste0("data/processed/", sp, "_model_data.rds"))
}

cat("\n=== DATA PREPARATION COMPLETE ===\n")
cat("\nCovariate structure:\n")
cat("  x.psi[i, 1:2]: intercept + water\n")
cat("  x.phi[i, k, 1:3]: intercept + clutter + water\n")
cat("  x.gamma[i, k, 1:3]: intercept + clutter + water\n")
cat("  x.p[i, k, j, 1:5]: intercept + clutter + temp + julian + moon\n")

# ===== DIAGNOSTIC CHECKS =====

cat("\n=== DIAGNOSTIC CHECKS ===\n")

test_data <- readRDS("data/processed/MYLU_model_data.rds")

cat("\nDimensions:\n")
cat("  y array:", paste(dim(test_data$y), collapse=" x "), "\n")
cat("  nsite:", test_data$nsite, "\n")
cat("  nyear:", test_data$nyear, "\n")

cat("\nExample y array (MYLU) - first 3 sites:\n")
for(i in 1:3) {
  cat("\nSite", i, "(", test_data$sites[i], "):\n")
  print(test_data$y[i, , ])
}

cat("\nVisits per site-year (J) distribution:\n")
print(table(test_data$J))

cat("\nNA check in covariates:\n")
cat("  x.psi NAs:", sum(is.na(test_data$x.psi)), "\n")
cat("  x.phi NAs:", sum(is.na(test_data$x.phi)), "\n")
cat("  x.gamma NAs:", sum(is.na(test_data$x.gamma)), "\n")
cat("  x.p NAs (surveyed occasions):", sum(is.na(test_data$x.p[!is.na(test_data$y)])), "\n")

