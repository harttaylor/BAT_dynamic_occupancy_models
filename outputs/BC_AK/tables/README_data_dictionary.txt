BAT OCCUPANCY MODEL OUTPUTS: ALASKA
Generated: 2026-09-29

The model was fitted to all BC and Alaska sites together, with region as a
fixed effect. These tables report the Alaska sites only. Covariate slopes
are shared by all sites; detection has its own Alaska offset. Each table has
an ALL_SPECIES_ version and per-species copies in per_species/. CIs are 95%.

site_year_occupancy    one row per species x site x year
  site              GRTS cell + quadrant, prefixed AK_
  surveyed          TRUE if the site was surveyed that year
  n_visits          visits (nights) used that year
  detected          1 = detected, 0 = not detected, NA = not surveyed
  p_occupied        posterior probability the site WAS occupied that year,
                    given what was observed there. 1 where detected, NA where
                    not surveyed. Use for 'where is it' maps. No CI: the state
                    is yes/no, so this probability is its full posterior.
  psi_expected      occupancy the model predicts from region, covariates and
                    dynamics, NOT conditioned on this site's own detections.
                    Exists for every site-year; use for 'how suitable' maps.

site_summary           one row per species x site
  mean_p_occupied   p_occupied averaged over surveyed years only
  mean_psi_expected psi_expected averaged over all years, with 95% CI

occupancy_trajectory   species x year: share of Alaska sites occupied
  naive             proportion of surveyed sites with a detection (raw data)
  n_surveyed        sites surveyed that year

dynamics               species x process x year (realised rates)
  Colonization      share of unoccupied sites in year t-1 occupied in year t
  Persistence       share of occupied sites in year t-1 still occupied in t
  Blank where too few sites were occupied to define the rate.

trends                 one row per species
  n_detections      detection nights: survey nights (visits) at a site with
                    at least one call identified to the species
  slope             linear change in occupancy per year
  change            last year minus first year
  p_increase        posterior probability last year > first year

covariate_effects      species x covariate. Logit scale, per 1 SD.
  Distances were log(x+1) transformed before standardising.
  Alaska vs BC      difference in per-visit detection at Alaska sites
  f                 posterior probability the effect has the sign of its mean
  excludes_zero     TRUE if the 95% CI does not include zero
