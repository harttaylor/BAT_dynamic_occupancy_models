################################################################################
# 04_audit_detections_by_region.R
#
# PURPOSE: Answer the question "are out-of-range occupancy estimates a model
# artifact, or is there something wrong with the data?"
#
# In this model y[i,t,j] ~ dbern(z[i,t] * p), so a single detection FORCES
# z = 1. There is no false-positive component. Therefore:
#
#   - Detections present in a region where the species doesn't occur
#       -> DATA problem (acoustic auto-ID false positives, or mislabeled sites)
#   - No detections, but occupancy still estimated
#       -> MODEL artifact (latent state filled in; gamma has no region term)
#
# Runs on the raw data only. No model re-run required.
################################################################################

library(dplyr)
library(tidyr)

RAW_DATA_FILE <- "data/raw/final_detection_data_combined.csv"
COVS_FILE     <- "data/processed/BCAK/site_covariates_final.csv"

START_YEAR <- 2017
END_YEAR   <- 2025
VALID_QUADRANTS <- c("NE", "NW", "SE", "SW")

SPECIES_CODES <- c("ANPA", "COTO", "EPFU", "EUMA", "LABO", "LACI", "LANO",
                   "MYCA", "MYCI", "MYEV", "MYLU", "MYSE", "MYTH", "MYVO",
                   "MYYU", "PAHE", "TABR")

# ============================================================================
# LOAD — mirror the filtering in 01c so site/region assignment matches
# ============================================================================

raw <- read.csv(RAW_DATA_FILE, stringsAsFactors = FALSE)
covs <- read.csv(COVS_FILE, stringsAsFactors = FALSE)

dat <- raw %>%
  mutate(
    site = paste0(GRTS.Cell.ID, "_", Quadrant),
    date = as.Date(Night),
    year = as.integer(format(date, "%Y"))
  ) %>%
  filter(Quadrant %in% VALID_QUADRANTS,
         year >= START_YEAR, year <= END_YEAR)

# Use the region from the COVARIATE file (this is what feeds region_idx in the
# model), not the Region column in the raw detection file. If these two
# disagree for any site, that is itself a finding — see check at bottom.
dat <- dat %>%
  inner_join(covs %>% dplyr::select(site, region_covfile = region), by = "site")

cat("Sites matched to covariate file:", n_distinct(dat$site), "\n\n")

# ============================================================================
# 1. DETECTION COUNTS BY SPECIES x REGION
# ============================================================================

det_long <- dat %>%
  dplyr::select(site, year, region = region_covfile, all_of(SPECIES_CODES)) %>%
  pivot_longer(all_of(SPECIES_CODES), names_to = "species", values_to = "count") %>%
  mutate(detected = !is.na(count) & count > 0)

# Nights with >=1 call sequence, and number of distinct sites involved
audit <- det_long %>%
  group_by(species, region) %>%
  summarize(
    nights_surveyed   = n(),
    nights_detected   = sum(detected),
    sites_surveyed    = n_distinct(site),
    sites_with_any_det = n_distinct(site[detected]),
    .groups = "drop"
  ) %>%
  mutate(pct_nights = round(100 * nights_detected / nights_surveyed, 2))

write.csv(audit, "outputs/audit_detections_by_species_region.csv", row.names = FALSE)

# Wide view: number of SITES with at least one detection, species x region
wide_sites <- audit %>%
  dplyr::select(species, region, sites_with_any_det) %>%
  pivot_wider(names_from = region, values_from = sites_with_any_det, values_fill = 0)

cat("=== SITES WITH >=1 DETECTION, by species x region ===\n")
print(as.data.frame(wide_sites))
cat("\n")

write.csv(wide_sites, "outputs/audit_sites_with_detections_wide.csv", row.names = FALSE)

# ============================================================================
# 2. FLAG SUSPECTED OUT-OF-RANGE DETECTIONS
# ============================================================================
# EDIT THIS LIST. These are the species/region combinations that should be
# biologically impossible. Add or remove based on your own range knowledge.

out_of_range <- tribble(
  ~species, ~region,
  "ANPA",   "Northern Region",
  "ANPA",   "Kootenay Region",
  "ANPA",   "South Coastal Region",
  "ANPA",   "Alaska",
  "EUMA",   "Northern Region",
  "EUMA",   "South Coastal Region",
  "EUMA",   "Alaska",
  "COTO",   "Northern Region",
  "COTO",   "Alaska",
  "TABR",   "Northern Region",
  "TABR",   "Alaska",
  "PAHE",   "Alaska",
  "MYSE",   "Alaska",
  "MYCI",   "Alaska",
  "MYTH",   "Alaska",
  "MYYU",   "Alaska",
  "LABO",   "Alaska"
)

flagged <- out_of_range %>%
  left_join(audit, by = c("species", "region")) %>%
  filter(!is.na(nights_surveyed))

cat("=== SUSPECTED OUT-OF-RANGE COMBINATIONS ===\n")
cat("Any row with nights_detected > 0 means the RAW DATA contains detections\n")
cat("where the species should not occur -> data / auto-ID problem.\n")
cat("Rows with nights_detected == 0 but occupancy in the model output\n")
cat("-> pure model artifact.\n\n")
print(as.data.frame(flagged))
cat("\n")

write.csv(flagged, "outputs/audit_out_of_range_flags.csv", row.names = FALSE)

# ============================================================================
# 3. THE OFFENDING RECORDS THEMSELVES
# ============================================================================
# If anything flagged above has detections, pull the actual rows so they can be
# checked against the recordings / auto-ID confidence.

offenders <- det_long %>%
  inner_join(out_of_range, by = c("species", "region")) %>%
  filter(detected) %>%
  arrange(species, region, site, year)

if (nrow(offenders) > 0) {
  cat("=== INDIVIDUAL OUT-OF-RANGE DETECTION RECORDS:", nrow(offenders), "===\n")
  print(head(as.data.frame(offenders), 50))
  write.csv(offenders, "outputs/audit_out_of_range_records.csv", row.names = FALSE)
  cat("\nFull list written to outputs/audit_out_of_range_records.csv\n\n")
} else {
  cat("No out-of-range detections in the raw data.\n")
  cat("=> Any out-of-range occupancy is a MODEL artifact, not a data problem.\n\n")
}

# ============================================================================
# 4. SENSITIVITY: does a stricter nightly threshold remove them?
# ============================================================================
# Current prep treats ONE call sequence as a detection. This is the most liberal
# possible threshold and the most exposed to auto-ID false positives.
# Check how much survives at >=2, >=3, >=5 sequences per night.

cat("=== EFFECT OF A STRICTER NIGHTLY DETECTION THRESHOLD ===\n")
for (thr in c(1, 2, 3, 5)) {
  n_flag <- det_long %>%
    inner_join(out_of_range, by = c("species", "region")) %>%
    filter(!is.na(count) & count >= thr) %>%
    nrow()
  cat(sprintf("  threshold >= %d call sequences/night: %d out-of-range detections remain\n",
              thr, n_flag))
}
cat("\nIf these drop sharply between 1 and 2-3, false-positive auto-ID is the\n")
cat("likely cause and a stricter threshold is a cheap, defensible fix.\n\n")

# ============================================================================
# 5. CONSISTENCY CHECK: region label in raw data vs covariate file
# ============================================================================

if ("Region" %in% names(dat)) {
  mismatch <- dat %>%
    distinct(site, Region, region_covfile) %>%
    filter(!is.na(Region), Region != "", Region != region_covfile)
  if (nrow(mismatch) > 0) {
    cat("WARNING:", nrow(mismatch), "sites have conflicting region labels\n")
    cat("between the raw detection file and the covariate file:\n")
    print(as.data.frame(mismatch))
    write.csv(mismatch, "outputs/audit_region_label_mismatches.csv", row.names = FALSE)
  } else {
    cat("Region labels consistent between raw data and covariate file.\n")
  }
}

# Duplicate site IDs across regions (possible GRTS cell ID collision BC vs AK)
dupes <- covs %>% count(site) %>% filter(n > 1)
if (nrow(dupes) > 0) {
  cat("\nWARNING: duplicate site IDs in covariate file:\n"); print(dupes)
} else {
  cat("No duplicate site IDs.\n")
}

cat("\nDone. Outputs in outputs/\n")
