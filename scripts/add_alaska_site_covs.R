library(dplyr)
library(elevatr)
library(sf)
library(terra)

# --- Load raw data, filter to AK ---
raw_data <- read.csv("data/raw/final_detection_data.csv")
sort(unique(raw_data$Region))  

ak_data <- raw_data %>%
  filter(Region == "Alaska")  

# --- Build one row per site ---
ak_sites <- ak_data %>%
  mutate(site = paste0(GRTS.Cell.ID, "_", Quadrant)) %>%
  filter(Quadrant %in% c("NE", "NW", "SE", "SW")) %>%
  group_by(site, GRTS.Cell.ID, Quadrant, Region) %>%
  summarize(
    lat  = mean(Lat,  na.rm = TRUE),   
    long = mean(Long, na.rm = TRUE),
    water_nearby = first(na.omit(Water.Nearby)),  # adjust if needed
    .groups = "drop"
  ) %>%
  rename(region = Region)

cat("AK sites:", nrow(ak_sites), "\n")


# --- Extract elevation using elevatr package ---
coords_sf <- st_as_sf(ak_sites, coords = c("long", "lat"), crs = 4326)
elev <- get_elev_point(coords_sf, src = "aws", z = 9)
ak_sites$elevation_m <- elev$elevation

# --- Distance to road (raster is already dist-to-road in meters) ---
roads <- rast("data/raw/roads_dist_AKonly.tif")  # already in EPSG:3338

ak_pts_3338 <- st_transform(coords_sf, crs = 3338)
ak_sites$dist_road_km <- extract(roads, vect(ak_pts_3338))[, 2] / 1000

# --- No harvest distance for AK — set NA, will be imputed to mean in prep script ---
ak_sites$dist_harvest_km <- NA

cat("Dist to road range:", round(min(ak_sites$dist_road_km, na.rm=TRUE), 1), "-",
    round(max(ak_sites$dist_road_km, na.rm=TRUE), 1), "km\n")

# --- Save ---
write.csv(ak_sites, "data/processed/ak_site_covariates.csv", row.names = FALSE)




# Look at distance to road file for AK 
# Read in file prepared by A. De Rosa 
roads <- rast("data/raw/roads_dist_AKonly.tif")
# Reproject to Alaska Albers for accurate meter-based distances
roads_proj <- project(roads, "EPSG:3338")

# Step 1: What does the raster actually look like?
roads_proj
res(roads_proj)
crs(roads_proj)
plot(roads_proj)  # visual check

# Step 2: What values are in it?
freq(roads_proj)        # all unique values + counts
minmax(roads_proj)      # min/max
hist(roads_proj)        # distribution

# Step 3: Do points fall within raster extent?
ext(roads_proj)   # raster bounding box
ak_pts_3338 <- st_transform(ak_pts, 3338)
st_bbox(ak_pts_3338)   # points bounding box — do these overlap?

# Quick sanity check — extents should overlap
ext(roads_proj)
st_bbox(ak_pts_3338)

# Extract distance values at each point
ak_sites$dist_road_m  <- extract(roads_proj, vect(ak_pts_3338))[, 2]
ak_sites$dist_road_km <- ak_sites$dist_road_m / 1000

# Check
summary(ak_sites$dist_road_km)
hist(ak_sites$dist_road_km, main = "AK sites: dist to road", xlab = "km")

# Which points are outside the raster extent?
pts_coords <- st_coordinates(ak_pts_3338)
rast_ext <- ext(roads_proj)

outside <- pts_coords[,1] < rast_ext[1] | pts_coords[,1] > rast_ext[2] |
  pts_coords[,2] < rast_ext[3] | pts_coords[,2] > rast_ext[4]

cat("Points outside raster extent:", sum(outside), "\n")
ak_sites[outside, c("site", "lat", "long")]

