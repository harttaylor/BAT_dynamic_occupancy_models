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
* builds and standardizes covariates (distances are log-transformed because of right-scewing)
* missing values are handled on a case-by-case basis
- *distance to harvest*: unsurveyed years filled with the site's nearest known year
- *clutter*: 40 sites had missing values in each year and were given the mean of all sites 
- *nightly mean temperature*: missing values filled from the site-year, then region-year, then region mean

### 01b_prepare_species_data.R
* builds each species detection history
* writes species_detection_summary.csv that is used in script 02 to fit models

### 02_fit_models.R
* fits one species at a time in top portion of script
* fits all species in a loop at bottom part 

### 03_figures_BC.R
* creates BC figures and tables 

### 03_figures_AK.R
* creates AK figures and tables

## Running the models 
* Each species can take hours to fit
* run_log_species.csv in the fits folder is updated after every species with run time and convergence (Rhat)

### Adding new data
1. Put the new data in data/raw/
2. Update the file names at the top of 00 and 01a, and END_YEAR in 01a
3. Run the scripts in order
4. Check any detections in a region where a species has not been recorded before: one detection is enough to include that region in the model

### Outputs 
* Iterations, run time and convergence for each species are recorded in run_log_species.csv in the fits folder
* *Occupancy, colonization, and persistence trends:* summaries are calculated from the estimated occupancy state of every site in every year
* Trends, each calculated in every posterior draw:
  - **Slope**: linear change in occupancy per year
  - **Change**: occupancy in 2025 minus 2017
  - **P(increase)**: the share of draws in 2025 where occupancy exceeds 2017. Trend figures colour species by this
* *Covariate effects* are on the logit scale per 1 SD of the covariate. For distances, a positive effect means the probability is higher farther from the feature
* *Site-level occupancy and maps:* maps colour each site by the model's estimate of whether it was occupied, which is corrected for imperfect detection. The symbol shape is raw data.
  - Fig 6 (site map): colour = estimated occupancy averaged over the years each site was surveyed; shape = ever detected
  - Fig 7 (map by year): surveyed site-years only; colour = that year's estimated occupancy; shape = detected that year.

