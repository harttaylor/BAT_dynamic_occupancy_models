BC BAT OCCUPANCY MODEL OUTPUTS -- DATA DICTIONARY
Generated: 2026-09-29

All files have one ALL_SPECIES_ version combining every species, and a
per-species version in per_species/.

Credible intervals are 95% (2.5% and 97.5% posterior quantiles) unless noted.

------------------------------------------------------------------------
ALL_SPECIES_site_year_occupancy.csv    one row per species x site x year
------------------------------------------------------------------------
  site              GRTS cell ID + quadrant, e.g. 110378_NE
  GRTS_Cell_ID      NABat GRTS cell
  Quadrant          NE / NW / SE / SW
  region            BC region
  lat, long         decimal degrees
  year
  surveyed          TRUE if the site was surveyed that year
  n_visits          number of visits (nights) used that year
  detected          1 = detected, 0 = not detected, NA = not surveyed
  p_occupied        posterior probability the site WAS occupied that year,
                    given what was observed there. Always 1 where detected.
                    NA where not surveyed. Use this for 'where is it' maps.
                    No CI: the underlying state is yes/no, so this single
                    probability is its full posterior.
  psi_expected      occupancy probability the model predicts for this site
                    from region, covariates and dynamics, NOT conditioned on
                    this site's own detections. Available for every
                    site-year. Use this for 'how suitable is it' maps or for
                    showing uncertainty.
  psi_expected_lo95, psi_expected_hi95   95% CI for psi_expected

------------------------------------------------------------------------
ALL_SPECIES_site_summary.csv    one row per species x site
------------------------------------------------------------------------
  n_years_surveyed, n_years_detected, ever_detected
  mean_p_occupied        p_occupied averaged over SURVEYED years only
  mean_psi_expected      psi_expected averaged over all years, with 95% CI

------------------------------------------------------------------------
ALL_SPECIES_occupancy_trajectory.csv    species x scope x year
------------------------------------------------------------------------
  scope   'BC (all sites)'        share of all sites occupied. Weighted by
                                  where sampling happened, so regions with
                                  more sites count for more.
          'BC (surveyed sites)'   same, restricted to sites surveyed that
                                  year. If this disagrees with 'all sites',
                                  the difference is model interpolation at
                                  unsurveyed site-years.
          'BC (region-weighted)'  average of the four regional values, so
                                  each region counts equally.
          <region name>           share of that region's sites occupied
  mean, lo95, hi95
  naive        proportion of surveyed sites with a detection (raw data)
  n_surveyed   sites surveyed that year in that scope

------------------------------------------------------------------------
ALL_SPECIES_dynamics.csv    species x scope x process x year
------------------------------------------------------------------------
  process   Colonization  share of sites unoccupied in year t-1 that were
                          occupied in year t
            Persistence   share of sites occupied in year t-1 still
                          occupied in year t
  Blank (NA) where too few sites were occupied to define the rate.

------------------------------------------------------------------------
ALL_SPECIES_trends.csv    species x scope
------------------------------------------------------------------------
  psi_first, psi_last    occupancy in first and last year
  slope                  linear change in occupancy per year, with 95% CI
  change                 last year minus first year, with 95% CI
  p_increase             posterior probability occupancy was higher in the
                         last year than the first
  n_detections           detection nights in that scope: survey nights (visits)
                         at a site with at least one call identified to the
                         species. Where 0, the species was never recorded
                         there: it is assumed absent, occupancy is fixed at
                         zero and there is no trend to report.

------------------------------------------------------------------------
ALL_SPECIES_covariate_effects.csv    species x covariate
------------------------------------------------------------------------
  Slopes on the logit scale, per 1 SD of the covariate. Shared across all
  regions (BC-wide). Distances were log(x+1) transformed before
  standardising.
  lo50, hi50       50% credible interval
  f                posterior probability the effect has the sign of its mean
  excludes_zero    TRUE if the 95% CI does not include zero
