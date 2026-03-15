################################################################################
# Add missing Alaska sites to the Occupancy tab export
#
#   1. Cami gave newly processed bat data (DataIDsOnly.xslx) from current year with all BC data
#      (old years + new year already combined)
#   2. That export was missing ~10 Alaska GRTS cells — this script adds them
#      back from BatActivity_2016to2024.csv
#   3. The output is the combined file used by all downstream scripts
#
# do NOT need to manually merge new-year BC data — that's already
# done in the Occupancy tab of the DataIDsOnly excel file. This script only fixes the missing Alaska sites.
################################################################################

library(dplyr)

# --- Input files (update occupancy_file path each year) ---
occupancy_file   <- "data/raw/final_detection_data.csv"    # latest Occupancy tab export
batactivity_file <- "data/raw/BatActivity_2016to2024.csv"  # original complete dataset

# --- Output ---
output_file <- "data/raw/final_detection_data_combined.csv"

# --- Load both ---
occ <- read.csv(occupancy_file, stringsAsFactors = FALSE)
bat <- read.csv(batactivity_file, stringsAsFactors = FALSE)

# --- Find Alaska GRTS cells missing from the Occupancy tab and add them ---
occ_ids     <- unique(occ$GRTS.Cell.ID)
bat_ak      <- bat %>% filter(Region == "Alaska")
missing_ids <- setdiff(unique(bat_ak$GRTS.Cell.ID), occ_ids)

ak_to_add <- bat_ak %>% filter(GRTS.Cell.ID %in% missing_ids)
combined  <- bind_rows(occ, ak_to_add)

cat(length(missing_ids), "Alaska GRTS cells added back\n")
cat("Combined:", nrow(combined), "rows,", n_distinct(combined$GRTS.Cell.ID), "GRTS cells\n")

write.csv(combined, output_file, row.names = FALSE)



################################################################################
# Extract Alaska site covariates
#
# Builds site covariates for Alaska from the combined data file.
# Extracts elevation (elevatr) and distance to road (raster).
# No distance-to-harvest for Alaska.
#
# Only needs re-running if new Alaska sites are added to the data.
#
# REQUIRES:
#   install.packages(c("elevatr", "sf", "terra"))
#   data/raw/roads_dist_AKonly.tif
#
# Output: data/ak_site_covariates.csv
################################################################################

library(dplyr)
library(elevatr)
library(sf)
library(terra)

# ============================================================================
# 1. BUILD ONE ROW PER ALASKA SITE-QUADRANT
# ============================================================================

raw_data <- read.csv("data/raw/final_detection_data_combined.csv")

ak_sites <- raw_data %>%
  filter(Region == "Alaska") %>%
  mutate(site = paste0(GRTS.Cell.ID, "_", Quadrant)) %>%
  filter(Quadrant %in% c("NE", "NW", "SE", "SW")) %>%
  group_by(site, GRTS.Cell.ID, Quadrant, Region) %>%
  summarize(
    lat  = mean(Lat,  na.rm = TRUE),
    long = mean(Long, na.rm = TRUE),
    water_nearby = first(na.omit(Water.Nearby)),
    .groups = "drop"
  ) %>%
  rename(grts_cell = GRTS.Cell.ID, region = Region)

cat("Alaska site-quadrants:", nrow(ak_sites), "\n")

# ============================================================================
# 2. EXTRACT ELEVATION
# ============================================================================

coords_sf <- st_as_sf(ak_sites, coords = c("long", "lat"), crs = 4326)
elev <- get_elev_point(coords_sf, src = "aws", z = 9)
ak_sites$elevation_m <- elev$elevation

cat("Elevation range:",
    round(min(ak_sites$elevation_m, na.rm = TRUE)), "-",
    round(max(ak_sites$elevation_m, na.rm = TRUE)), "m\n")

# ============================================================================
# 3. EXTRACT DISTANCE TO ROAD
# ============================================================================

roads <- rast("data/raw/roads_dist_AKonly.tif")
ak_pts_3338 <- st_transform(coords_sf, crs = 3338)
ak_sites$dist_road_km <- extract(roads, vect(ak_pts_3338))[, 2] / 1000

cat("Dist to road range:",
    round(min(ak_sites$dist_road_km, na.rm = TRUE), 1), "-",
    round(max(ak_sites$dist_road_km, na.rm = TRUE), 1), "km\n")

# No harvest distance for Alaska
ak_sites$dist_harvest_km <- NA

# ============================================================================
# 4. SAVE
# ============================================================================

write.csv(ak_sites, "data/ak_site_covariates.csv", row.names = FALSE)
cat("Saved:", nrow(ak_sites), "Alaska sites to data/ak_site_covariates.csv\n")
