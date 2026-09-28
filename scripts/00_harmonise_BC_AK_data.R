################################################################################
# 00_harmonise_BC_AK.R
#
# ALASKA REPORT ONLY. Skip this script for the BC report.
#
# Puts the Alaska file into the same format as the BC file and stacks the two
# into one file, so the rest of the scripts can treat Alaska as a fifth region.
#
# OUTPUT: data/processed/AllYearActivitybyNight_BC_AK_combined.csv
################################################################################

library(dplyr)

# ============================================================================
# SETTINGS  (update the file names when new data arrive)
# ============================================================================

BC_FILE  <- "data/raw/AllYearActivitybyNight_2016_2025_BC_withWeather_NEAR.csv"
AK_FILE  <- "data/raw/All Year Activity by Night 2016-2025 AK-Aug19.csv"
OUT_FILE <- "data/processed/AllYearActivitybyNight_BC_AK_combined.csv"

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# 1. LOAD AND CLEAN COLUMN NAMES
# ============================================================================
# Spaces, dots and brackets in column names all become "_", in both files.

bc <- read.csv(BC_FILE)
ak <- read.csv(AK_FILE)

names(bc) <- gsub("_+$", "", gsub("[^A-Za-z0-9]+", "_", names(bc)))
names(ak) <- gsub("_+$", "", gsub("[^A-Za-z0-9]+", "_", names(ak)))

# The AK file is the source for Alaska, so drop any Alaska rows in the BC file
bc <- bc %>% filter(Region != "Alaska")

# ============================================================================
# 2. MAKE THE AK FILE MATCH THE BC FILE
# ============================================================================

# Rename AK columns to the BC names. Written as  BC name = "AK name".
# any_of() means this line can safely be run twice.
ak <- ak %>%
  rename(any_of(c(Latitude       = "Lat",
                  Longitude      = "Long",
                  DIST_ROAD      = "road_dist",
                  DIST_HARVEST   = "harv_dist",
                  YEAR_HARVEST   = "harv_year",
                  Flyingsquirrel = "lyingsquirrel")))

# Frequency columns: R reads "25K" as "X25K"; BC calls it "F25K"
names(ak) <- sub("^X([0-9]+K)$", "F\\1", names(ak))

# Distance to water: BC has one distance to the nearest water; AK has streams
# and lakes separately, so use whichever is nearer
ak$DIST_WATER_M <- pmin(ak$strm_dist, ak$lake_dist, na.rm = TRUE)

ak$Region <- "Alaska"

# Same date format in both files (YYYY-MM-DD)
bc$Night <- format(as.Date(bc$Night), "%Y-%m-%d")
ak$Night <- format(as.Date(ak$Night), "%Y-%m-%d")

# ============================================================================
# 3. CHECKS: both of these should print character(0), i.e. nothing
# ============================================================================

# AK species columns that BC doesn't have (they would be lost when stacking)
setdiff(grep("^[A-Z]{4,}$", names(ak), value = TRUE), names(bc))

# Columns the models need that AK is missing
setdiff(c("GRTS_Cell_ID", "Quadrant", "Night", "Region", "Latitude", "Longitude",
          "Elevation", "Percent_Clutter", "DIST_WATER_M", "DIST_ROAD",
          "DIST_HARVEST", "YEAR_HARVEST", "Nightly_Mean_Temp"), names(ak))

# ============================================================================
# 4. STACK AND SAVE
# ============================================================================
# Only the BC columns are kept. Species AK never records (e.g. EPFU) are left
# empty (NA) for Alaska rows; 01b treats that as "not detected".

ak <- ak[, intersect(names(bc), names(ak))]

combined <- bind_rows(bc, ak)

table(combined$Region)   # rows per region

write.csv(combined, OUT_FILE, row.names = FALSE)
