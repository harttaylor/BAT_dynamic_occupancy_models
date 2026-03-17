################################################################################
# Per-Species Regional Occupancy with CIs
#
# For each species, creates a faceted figure (one panel per region) showing
# model-estimated occupancy with 95% credible intervals.
#
# CIs are computed from the full posterior of z (fit$sims.list$z).
# For each MCMC iteration, we compute the mean z across sites in a region,
# then take quantiles across iterations.
#
# REQUIRES: z in parameters.to.save when fitting models (02a/02b).
#
# Output: outputs/[MODE]/figures/[SPECIES]_04_regional.png
################################################################################

library(ggplot2)
library(dplyr)

# ============================================================================
# SETTINGS
# ============================================================================

MODE <- "BC"   # "BC" or "BCAK"

FITS_DIR <- file.path("outputs", MODE, "fits")
DATA_DIR <- file.path("data/processed", MODE)
FIG_DIR  <- file.path("outputs", MODE, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

MODE_LABEL <- ifelse(MODE == "BC", "BC", "BC + Alaska")

region_colors <- c(
  "Kootenay Region"      = "#1B9E77",
  "Northern Region"      = "#D95F02",
  "Okanagan Region"      = "#7570B3",
  "South Coastal Region" = "#E7298A",
  "Alaska"               = "#66A61E"
)

# Fill colors (lighter versions for CI ribbons)
region_fills <- c(
  "Kootenay Region"      = "#1B9E77",
  "Northern Region"      = "#D95F02",
  "Okanagan Region"      = "#7570B3",
  "South Coastal Region" = "#E7298A",
  "Alaska"               = "#66A61E"
)

# ============================================================================
# FIND SPECIES
# ============================================================================

fit_files <- list.files(FITS_DIR, pattern = "QUAD_.*_fit\\.rds")
species_list <- gsub("QUAD_|_fit\\.rds", "", fit_files)
cat("Species:", paste(species_list, collapse = ", "), "\n\n")

# ============================================================================
# LOOP THROUGH SPECIES
# ============================================================================

for (sp in species_list) {
  
  cat(sp, "... ")
  
  fit <- readRDS(file.path(FITS_DIR, paste0("QUAD_", sp, "_fit.rds")))
  d   <- readRDS(file.path(DATA_DIR, paste0("QUAD_", sp, "_data.rds")))
  
  years        <- d$years
  nyear        <- d$nyear
  region_names <- d$regions
  region_idx   <- d$region
  
  # fit$sims.list$z is an array [n_iter x nsite x nyear]
  z_post <- fit$sims.list$z
  n_iter <- dim(z_post)[1]
  
  # Build regional occupancy data frame with CIs
  reg_occ <- data.frame()
  
  for (r in seq_along(region_names)) {
    sites_r <- which(region_idx == r)
    n_r     <- length(sites_r)
    if (n_r == 0) next
    
    for (t in 1:nyear) {
      # Extract z values for this region and year: matrix [n_iter x n_sites_in_region]
      if (n_r == 1) {
        # Single site: z_post[, sites_r, t] is a vector of length n_iter
        iter_means <- z_post[, sites_r, t]
      } else {
        # Multiple sites: z_post[, sites_r, t] is a matrix [n_iter x n_sites]
        # Take the mean across sites for each iteration
        iter_means <- rowMeans(z_post[, sites_r, t])
      }
      
      reg_occ <- rbind(reg_occ, data.frame(
        region_label = paste0(region_names[r], " (n = ", n_r, ")"),
        region       = region_names[r],
        year         = years[t],
        psi          = mean(iter_means),
        lcl          = quantile(iter_means, 0.025),
        ucl          = quantile(iter_means, 0.975),
        n_sites      = n_r
      ))
    }
  }
  
  yr_label <- paste0(min(years), "\u2013", max(years))
  
  p <- ggplot(reg_occ, aes(x = year, y = psi, color = region, fill = region)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    facet_wrap(~ region_label, ncol = 2) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_x_continuous(breaks = years) +
    scale_color_manual(values = region_colors, guide = "none") +
    scale_fill_manual(values = region_fills, guide = "none") +
    labs(
      title = paste0(sp, " \u2014 Regional Occupancy (", MODE_LABEL, ", ", yr_label, ")"),
      subtitle = "Model-estimated occupancy per region | 95% credible intervals",
      x = "Year", y = "Occupancy (\u03c8)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "gray30"),
      strip.text = element_text(face = "bold", size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
    )
  
  ggsave(file.path(FIG_DIR, paste0(sp, "_04_regional.png")),
         p, width = 12, height = 8, dpi = 300)
  
  cat("done\n")
}

cat("\nRegional figures saved to:", FIG_DIR, "\n")
