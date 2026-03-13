################################################################################
# BC Bat Monitoring - Add Elevation from DEM using elevatr
#
# RUN THIS AFTER data prep to add real elevation values
# Then re-run data prep with elevation included
#
# REQUIRES: 
#   install.packages("elevatr")
#   install.packages("sf")
#
# This script:
# 1. Loads site coordinates
# 2. Fetches elevation from AWS terrain tiles (SRTM-based)
# 3. Saves updated site info with elevation
#
################################################################################

library(dplyr)
library(elevatr)
library(sf)

# ============================================================================
# SETTINGS
# ============================================================================

INPUT_FILE <- "data/processed/bc_only_2017start/site_covariates_v5.csv"
OUTPUT_FILE <- "data/processed/bc_only_2017start/site_covariates_with_elevation.csv"

# ============================================================================
# LOAD SITE DATA
# ============================================================================

site_data <- read.csv(INPUT_FILE)

cat("Loaded", nrow(site_data), "sites\n")
cat("Sites with lat/long:", sum(!is.na(site_data$lat) & !is.na(site_data$long)), "\n\n")

# Check for missing coordinates
missing_coords <- is.na(site_data$lat) | is.na(site_data$long)
if (any(missing_coords)) {
  cat("WARNING:", sum(missing_coords), "sites missing coordinates\n")
  cat("These will get NA elevation\n\n")
}

# ============================================================================
# FETCH ELEVATION
# ============================================================================

cat("Fetching elevation from AWS terrain tiles...\n")
cat("(This may take a minute for first download)\n\n")

# Create spatial points
coords_df <- site_data %>%
  filter(!is.na(lat), !is.na(long)) %>%
  select(site, long, lat)

# Convert to sf object (required by elevatr)
coords_sf <- st_as_sf(coords_df, coords = c("long", "lat"), crs = 4326)

# Get elevation - uses AWS terrain tiles (SRTM-based, ~30m resolution)
# z parameter controls zoom level: higher = more detail but slower
# z=9 is usually good balance for provincial-scale data
elevation_data <- get_elev_point(coords_sf, src = "aws", z = 9)

# Extract elevation values
coords_df$elevation_m <- elevation_data$elevation

cat("Elevation retrieved for", sum(!is.na(coords_df$elevation_m)), "sites\n")
cat("Elevation range:", round(min(coords_df$elevation_m, na.rm = TRUE)), "-", 
    round(max(coords_df$elevation_m, na.rm = TRUE)), "m\n\n")

# ============================================================================
# MERGE BACK TO SITE DATA
# ============================================================================

site_data <- site_data %>%
  left_join(coords_df %>% select(site, elevation_m), by = "site")

# Summary
cat("Elevation summary:\n")
print(summary(site_data$elevation_m))

# ============================================================================
# SAVE
# ============================================================================

write.csv(site_data, OUTPUT_FILE, row.names = FALSE)



