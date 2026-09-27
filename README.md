# Dynamic Occupancy Analysis of Bat Acoustic Detection Data 
Multi-season dynamic occupancy analysis using 9 years of bat acoustic detection data across British Columbia (BC) and Alaska (AK). Models estimate species-level trends in occupancy, colonization, and persistence while accounting for imperfect detection. Includes evaluation of optimal sampling deisgn scenarios to inform future monitoring efforts.

## Scripts 
### 00_harmonise_BC_AK.R
* renames the AK columns to the BC names, e.g., Lat -> Latitude, harv_dist → DIST_HARVEST, road_dist → DIST_ROAD, X40K → F40K.
* builds DIST_WATER_M for AK as the nearer of stream and lake distance, which matches AK's own WaterDist
* writes one combined file with the BC column layout (OUTPUT: data/processed/AllYearActivitybyNight_BC_AK_combined.csv)