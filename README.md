# Bat population dynamics and trends across British Columbia and Alaska
Bayesian dynamic occupancy analysis via JAGS of 9 years (2017-2025) of acoustic bat detection data across British Columbia (BC) + Alaska (AK). Models estimate species-level trends in occupancy, colonization, and persistence while accounting for imperfect detection. Includes evaluation of optimal sampling design scenarios to inform future monitoring efforts.

Two separate analyses were run: 1) BC data only, and 2) BC + AK data. The BC analysis models 185 BC sites in four regions and 15 bat species. The Alaska analysis adds 15 Alaska sites as a fifth region and models the 7 species detected in Alaska. AK is modeled with BC because 15 sites are too few to fit the model alone - differences in model structure and procedure for running BC vs. BC + AK are below. 

### Data and study design
The analysis uses nightly acoustic detection data from passive detectors at NABat GRTS cells in BC, 2017-2025. 2016 was excluded because few sites were surveyed.
* *Site:* one GRTS cell quadrant. Only the four standard quadrants (NE, NW, SE, SW) were included.
* *Visit:* one survey night. Deployments run 7-8 consecutive nights; we retained every second night to reduce dependence between consecutive nights, with max 4 visits per site-year.
A site year needed minimum 2 visits (to inform detection), and a site needed at least 2 surveyed years (to inform colonization and persistence).



### Model overview 
Each species was fitted separately with a Bayesian dynamic occupancy model. Initial occupancy was modelled as a function of distance to water and elevation, colonization and persistence both include distance to water, elevation, year-varying distance to harvest, distance to road (does not vary yearly).

Detection has percentage clutter, temperature, and Julian date (all vary at the site-night level), and a random year effect to absorb changes in classifier, hardware, vetting over the 9-year period, and prevent confounding linear trends in the detection process (p) with occupancy trends.

*Region* is a fixed effect on initial occupancy, colonization, and persistence, with separate intercepts for each region. 

*Year effects* are random (normal, mean 0, estimated SD) on colonization, persistence, and detection, and are shared across regions. The detection year effect absorbs year-to-year changes in detectability (classifiers, hardware, vetting) so they are not interpreted as changes in occupancy.

*Accounting for varying species ranges:* For each species, occupancy is fixed at zero in any region where it had no detection nights in 2017-2025, equivalent to restricting that species' analysis to regions where it was recorded. 
* Without this constraint the model has no data to rule out "occupied in 2017, then lost" in those regions, resulting in 2017 occupancy being set by the prior rather than the data, producing spurious declines and distorted BC-wide trends. 
* This constraint on occupancy works in conjunction with the region effects on colonization and persistence and is automatically applied each time the models are run, so a species first detected in a new region in a future year will be modeled there, *allowing the model to capture potential new range expansions*. 

### To run the models with BC data only skip script 00, and just make sure to set the input and output to BC folders at the top of scripts 01a-02, and create figures using 03_figures_BC

### To run the models for Alaska + BC combined start with 00, and set input/output directories to AK at the top of scripts 01a-02, and run 03_figures_AK to make figures using the Alaska region only
The differences between these scripts are:
1. When modelling Alaska we add an extra level to the region fixed effect, and Alaska has its own intercepts for initial occupancy, colonization and persistence. 
2. Only species that are detected in Alaska are modeled and thus have resulting outputs and figures made. This is done automatically and based on detection histories in 01b_prepare_species_data, so if a species range ever expands to Alaska it will be detected and modeled. 
3. *Alaska detection offset:* per-visit detection at Alaska sites includes a fixed offset on the logit scale, so detectability can differ between Alaska and BC (accounting for different recording equipment, auto-ID classification, shorter summer nights etc.)


## Model requirements
* R packages: dplyr, tidyr, ggplot2, jagsUI
* JAGS installed separately: https://sourceforge.net/projects/mcmc-jags/


## Scripts 
### 00_harmonise_BC_AK.R
* renames the AK columns to the BC names, e.g., Lat -> Latitude, harv_dist → DIST_HARVEST, road_dist → DIST_ROAD, X40K → F40K.
* builds DIST_WATER_M for AK as the nearer of stream and lake distance, which matches AK's own WaterDist
* writes one combined file with the BC column layout (OUTPUT: data/processed/AllYearActivitybyNight_BC_AK_combined.csv)

### 01a_prepare_covariates.R
* applies the data rules
* builds and standardises covariates 
