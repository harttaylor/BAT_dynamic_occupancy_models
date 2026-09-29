# Bat population dynamics and trends across British Columbia and Alaska
Bayesian dynamic occupancy analysis via JAGS of 9 years (2017-2025) of acoustic bat detection data across British Columbia (BC) + Alaska (AK). Models estimate species-level trends in occupancy, colonization, and persistence while accounting for imperfect detection. Includes evaluation of optimal sampling design scenarios to inform future monitoring efforts.

### Model overview 
Models initial occupancy as a function of distance to water and elevation, colonization and persistence both include distance to water, elevation, year-varying distance to harvest, distance to road (does not vary yearly).

Detection has percentage clutter, temperature, and Julian date (all vary at the site-night level), and a random year effect to absorb changes in classifier, hardware, vetting over the 9-year period, and prevent confounding linear trends in the detection process (p) with occupancy trends.

Region is a fixed effect on initial occupancy, colonization, and persistence, and occupancy is initialized according to whether the species was ever detected in a given region.
- A species is assumed ABSENT from any region where it was never detected (present[r] = 0): its occupancy there is fixed at zero. 
- Without this, the model has no data to rule out "occupied in 2017, then lost", so the first
year in those regions is set by the prior and shows a spurious decline. 
This allows us to 
- *a)* make region-level figures and show differences in occupancy, colonization, and persistence trends across different regions 
- *b)* avoid having the same occupancy starting place across regions (initial occupancy prior uniform across regions), making it look like occupancy was high in year 1 and then lost (showing a spurious decline)

### To run the models with BC data only skip script 00, and just make sure to set the input and output to BC folders at the top of scripts 01a-02, and create figures using 03_figures_BC

### To run the models for Alaska + BC combined start with 00, and set input/output directories to AK at the top of scripts 01a-02, and run 03_figures_AK to make figures using the Alaska region only
The differences between these scripts are:
1. When modelling Alaska we add an extra level to the region fixed effect
2. Only species that are detected in Alaska are modelled and thus have resulting outputs and figures made. This is done automatically and based on detection histories in 01b_prepare_species_data, so if a species range ever expands to Alaska it will be detected and modelled. 

## Scripts 
### 00_harmonise_BC_AK.R
* renames the AK columns to the BC names, e.g., Lat -> Latitude, harv_dist → DIST_HARVEST, road_dist → DIST_ROAD, X40K → F40K.
* builds DIST_WATER_M for AK as the nearer of stream and lake distance, which matches AK's own WaterDist
* writes one combined file with the BC column layout (OUTPUT: data/processed/AllYearActivitybyNight_BC_AK_combined.csv)

