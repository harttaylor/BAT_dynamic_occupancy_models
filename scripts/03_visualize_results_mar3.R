################################################################################
# 03_visualize_results.R — BC Bat Dynamic Occupancy Figures + Text Summaries
#
# Creates figures and brief written results from fitted models (02_fit_models.R).
#
# Set SPECIES below to make outputs for one species, or set RUN_ALL <- TRUE
# to loop through all available fits.
#
# Per-species outputs:
#   1. Occupancy trend plot
#   2. Year-specific colonization & persistence plot
#   3. Covariate effects forest plot
#   4. Brief text summary (results.txt)
#
# Multi-species outputs (when RUN_ALL = TRUE):
#   5. All-species occupancy panel
#   6. All-species λ comparison
#   7. All-species covariate effects
#   8. All-species dynamics panel
#   9. Summary CSV
################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)

# ============================================================================
# SETTINGS
# ============================================================================

SPECIES  <- "TABR"     # Single species (ignored if RUN_ALL = TRUE)
RUN_ALL  <- TRUE     # Set TRUE to process all available species

INPUT_DIR  <- "outputs/bc_dynocc"
DATA_DIR   <- "data/processed"

FIG_DIR    <- "outputs/bc_dynocc/figures"

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# HELPER: EXTRACT RESULTS FROM A FIT
# ============================================================================

extract_results <- function(species, input_dir, data_dir) {
  
  fit_file <- file.path(input_dir, paste0("QUAD_", species, "_fit.rds"))
  if (!file.exists(fit_file)) return(NULL)
  
  fit <- readRDS(fit_file)
  summ <- fit$summary
  
  # Load data for metadata
  d <- readRDS(file.path(data_dir, paste0("QUAD_", species, "_data.rds")))
  years <- d$years
  nyear <- d$nyear
  
  # Annual occupancy
  occ <- data.frame(
    species = species,
    year = years,
    psi = summ[paste0("psi.fs[", 1:nyear, "]"), "mean"],
    lcl = summ[paste0("psi.fs[", 1:nyear, "]"), "2.5%"],
    ucl = summ[paste0("psi.fs[", 1:nyear, "]"), "97.5%"]
  )
  
  # Year-specific dynamics (nyear-1 transitions)
  # X-axis = transition endpoint year (year 2 through final year)
  dyn <- data.frame(
    species = species,
    year_from = years[1:(nyear-1)],
    year_to   = years[2:nyear],
    gamma_mean = summ[paste0("gamma.avg[", 1:(nyear-1), "]"), "mean"],
    gamma_lcl  = summ[paste0("gamma.avg[", 1:(nyear-1), "]"), "2.5%"],
    gamma_ucl  = summ[paste0("gamma.avg[", 1:(nyear-1), "]"), "97.5%"],
    phi_mean   = summ[paste0("phi.avg[", 1:(nyear-1), "]"), "mean"],
    phi_lcl    = summ[paste0("phi.avg[", 1:(nyear-1), "]"), "2.5%"],
    phi_ucl    = summ[paste0("phi.avg[", 1:(nyear-1), "]"), "97.5%"],
    ext_mean   = summ[paste0("ext.avg[", 1:(nyear-1), "]"), "mean"],
    ext_lcl    = summ[paste0("ext.avg[", 1:(nyear-1), "]"), "2.5%"],
    ext_ucl    = summ[paste0("ext.avg[", 1:(nyear-1), "]"), "97.5%"]
  )
  
  # Trend summary
  trend <- data.frame(
    species = species,
    nsite = d$nsite,
    nyear = nyear,
    year_range = paste0(min(years), "-", max(years)),
    lambda_mean = summ["lambda", "mean"],
    lambda_lcl  = summ["lambda", "2.5%"],
    lambda_ucl  = summ["lambda", "97.5%"],
    slope_mean  = summ["trend.slope", "mean"],
    slope_lcl   = summ["trend.slope", "2.5%"],
    slope_ucl   = summ["trend.slope", "97.5%"],
    mean_gamma  = summ["mean.gamma", "mean"],
    mean_gamma_lcl = summ["mean.gamma", "2.5%"],
    mean_gamma_ucl = summ["mean.gamma", "97.5%"],
    mean_phi    = summ["mean.phi", "mean"],
    mean_phi_lcl = summ["mean.phi", "2.5%"],
    mean_phi_ucl = summ["mean.phi", "97.5%"],
    mean_ext    = summ["mean.ext", "mean"],
    mean_p      = summ["mean.p", "mean"],
    mean_p_lcl  = summ["mean.p", "2.5%"],
    mean_p_ucl  = summ["mean.p", "97.5%"],
    max_rhat    = max(summ[, "Rhat"], na.rm = TRUE),
    min_neff    = min(summ[, "n.eff"], na.rm = TRUE)
  )
  
  # Covariate effects (non-intercept)
  coef <- data.frame(
    species = species,
    parameter = c(
      "ψ: Elevation", "ψ: Water Nearby",
      "γ: Dist to Road", "γ: Dist to Harvest",
      "φ: Dist to Road", "φ: Dist to Harvest",
      "p: Clutter", "p: Temperature", "p: Julian"
    ),
    category = c(
      rep("Initial Occupancy", 2),
      rep("Colonization", 2),
      rep("Persistence", 2),
      rep("Detection", 3)
    ),
    estimate = c(
      summ["beta.psi[2]", "mean"], summ["beta.psi[3]", "mean"],
      summ["beta.gamma[1]", "mean"], summ["beta.gamma[2]", "mean"],
      summ["beta.phi[1]", "mean"], summ["beta.phi[2]", "mean"],
      summ["beta.p[2]", "mean"], summ["beta.p[3]", "mean"], summ["beta.p[4]", "mean"]
    ),
    lcl = c(
      summ["beta.psi[2]", "2.5%"], summ["beta.psi[3]", "2.5%"],
      summ["beta.gamma[1]", "2.5%"], summ["beta.gamma[2]", "2.5%"],
      summ["beta.phi[1]", "2.5%"], summ["beta.phi[2]", "2.5%"],
      summ["beta.p[2]", "2.5%"], summ["beta.p[3]", "2.5%"], summ["beta.p[4]", "2.5%"]
    ),
    ucl = c(
      summ["beta.psi[2]", "97.5%"], summ["beta.psi[3]", "97.5%"],
      summ["beta.gamma[1]", "97.5%"], summ["beta.gamma[2]", "97.5%"],
      summ["beta.phi[1]", "97.5%"], summ["beta.phi[2]", "97.5%"],
      summ["beta.p[2]", "97.5%"], summ["beta.p[3]", "97.5%"], summ["beta.p[4]", "97.5%"]
    )
  )
  coef$significant <- ifelse(coef$lcl > 0 | coef$ucl < 0, "Yes", "No")
  
  return(list(occ = occ, dyn = dyn, trend = trend, coef = coef, years = years))
}

# ============================================================================
# COLOR PALETTE
# ============================================================================

cat_colors <- c(
  "Initial Occupancy" = "#1B9E77",
  "Colonization" = "#D95F02",
  "Persistence" = "#7570B3",
  "Detection" = "#666666"
)

# ============================================================================
# DETERMINE SPECIES TO PROCESS
# ============================================================================

if (RUN_ALL) {
  fit_files <- list.files(INPUT_DIR, pattern = "QUAD_.*_fit\\.rds", full.names = FALSE)
  species_list <- gsub("QUAD_|_fit\\.rds", "", fit_files)
  cat("Processing all species:", paste(species_list, collapse = ", "), "\n\n")
} else {
  species_list <- SPECIES
}

# ============================================================================
# EXTRACT ALL RESULTS
# ============================================================================

all_occ   <- list()

all_dyn   <- list()
all_trend <- list()
all_coef  <- list()

for (sp in species_list) {
  cat("Loading", sp, "...")
  res <- extract_results(sp, INPUT_DIR, DATA_DIR)
  if (is.null(res)) { cat(" NOT FOUND\n"); next }
  all_occ[[sp]]   <- res$occ
  all_dyn[[sp]]   <- res$dyn
  all_trend[[sp]] <- res$trend
  all_coef[[sp]]  <- res$coef
  cat(" OK (λ =", round(res$trend$lambda_mean, 2), ")\n")
}

occ_df   <- bind_rows(all_occ)


dyn_df   <- bind_rows(all_dyn)
trend_df <- bind_rows(all_trend)
coef_df  <- bind_rows(all_coef)

cat("\nLoaded", length(all_occ), "species\n\n")

# ============================================================================
# PER-SPECIES FIGURES + TEXT SUMMARIES
# ============================================================================

cat("Creating figures...\n")

param_order <- c(
  "ψ: Elevation", "ψ: Water Nearby",
  "γ: Dist to Road", "γ: Dist to Harvest",
  "φ: Dist to Road", "φ: Dist to Harvest",
  "p: Clutter", "p: Temperature", "p: Julian"
)

for (sp in species_list) {
  
  sp_occ   <- occ_df   %>% filter(species == sp)
  sp_dyn   <- dyn_df   %>% filter(species == sp)
  sp_trend <- trend_df %>% filter(species == sp)
  sp_coef  <- coef_df  %>% filter(species == sp)
  
  year_range <- paste0(min(sp_occ$year), "\u2013", max(sp_occ$year))
  
  # --- Figure 1: Occupancy trend ---
  p1 <- ggplot(sp_occ, aes(x = year, y = psi)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.25, fill = "#1B9E77") +
    geom_line(linewidth = 1.3, color = "#1B9E77") +
    geom_point(size = 3.5, color = "#1B9E77") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_x_continuous(breaks = sp_occ$year) +
    labs(
      title = paste0(sp, " \u2014 Occupancy Trend (BC, ", year_range, ")"),
      subtitle = sprintf("\u03bb = %.2f (%.2f, %.2f) | Mean p = %.2f",
                          sp_trend$lambda_mean, sp_trend$lambda_lcl,
                          sp_trend$lambda_ucl, sp_trend$mean_p),
      x = "Year", y = "Occupancy (\u03c8)"
    ) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(FIG_DIR, paste0(sp, "_01_occupancy_trend.png")),
         p1, width = 10, height = 6, dpi = 300)
  
  # --- Figure 2: Year-specific dynamics ---
  # X-axis shows all transition endpoint years
  dyn_years <- sp_dyn$year_to
  
  p2 <- ggplot(sp_dyn) +
    geom_ribbon(aes(x = year_to, ymin = gamma_lcl, ymax = gamma_ucl),
                alpha = 0.2, fill = "#D95F02") +
    geom_line(aes(x = year_to, y = gamma_mean, color = "Colonization"), linewidth = 1.2) +
    geom_point(aes(x = year_to, y = gamma_mean, color = "Colonization"), size = 3) +
    geom_ribbon(aes(x = year_to, ymin = phi_lcl, ymax = phi_ucl),
                alpha = 0.2, fill = "#7570B3") +
    geom_line(aes(x = year_to, y = phi_mean, color = "Persistence"), linewidth = 1.2) +
    geom_point(aes(x = year_to, y = phi_mean, color = "Persistence"), size = 3) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_x_continuous(breaks = dyn_years) +
    scale_color_manual(values = c("Colonization" = "#D95F02", "Persistence" = "#7570B3")) +
    labs(
      title = paste0(sp, " \u2014 Year-Specific Dynamics (BC, ",
                     min(dyn_years), "\u2013", max(dyn_years), ")"),
      subtitle = "At average covariate values | 95% credible intervals",
      x = "Year (transition endpoint)", y = "Probability", color = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
  
  ggsave(file.path(FIG_DIR, paste0(sp, "_02_dynamics.png")),
         p2, width = 10, height = 6, dpi = 300)
  
  # --- Figure 3: Covariate effects ---
  sp_coef$parameter <- factor(sp_coef$parameter, levels = rev(param_order))
  
  p3 <- ggplot(sp_coef, aes(x = estimate, y = parameter, color = category)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = lcl, xmax = ucl), height = 0.35, linewidth = 0.9) +
    geom_point(aes(shape = significant), size = 3.5) +
    scale_color_manual(values = cat_colors, name = "Process") +
    scale_shape_manual(values = c("Yes" = 16, "No" = 1), name = "95% CI\nexcludes 0") +
    labs(
      title = paste0(sp, " \u2014 Covariate Effects (BC, ", year_range, ")"),
      subtitle = "Standardized effects on logit scale | 95% credible intervals",
      x = "Effect size", y = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "right", plot.title = element_text(face = "bold"))
  
  ggsave(file.path(FIG_DIR, paste0(sp, "_03_covariate_effects.png")),
         p3, width = 10, height = 7, dpi = 300)
  
  # --- Text summary ---
  sig_covs <- sp_coef %>% filter(significant == "Yes")
  
  trend_word <- if (sp_trend$lambda_ucl < 1) {
    "DECLINE"
  } else if (sp_trend$lambda_lcl > 1) {
    "INCREASE"
  } else {
    "STABLE (no clear trend)"
  }
  
  txt <- c(
    paste0("RESULTS: ", sp),
    paste0(rep("=", 50)),
    "",
    paste0("Data: ", sp_trend$nsite, " sites, ", sp_trend$nyear, " years (",
           sp_trend$year_range, ")"),
    "",
    "TREND",
    sprintf("  lambda = %.2f (%.2f, %.2f) -> %s",
            sp_trend$lambda_mean, sp_trend$lambda_lcl, sp_trend$lambda_ucl, trend_word),
    sprintf("  Trend slope = %.4f (%.4f, %.4f)",
            sp_trend$slope_mean, sp_trend$slope_lcl, sp_trend$slope_ucl),
    sprintf("  Occupancy: %.2f (%d) -> %.2f (%d)",
            sp_occ$psi[1], sp_occ$year[1],
            sp_occ$psi[nrow(sp_occ)], sp_occ$year[nrow(sp_occ)]),
    "",
    "DYNAMICS (mean across years, at avg covariate values)",
    sprintf("  Colonization: %.3f (%.3f, %.3f)",
            sp_trend$mean_gamma, sp_trend$mean_gamma_lcl, sp_trend$mean_gamma_ucl),
    sprintf("  Persistence:  %.3f (%.3f, %.3f)",
            sp_trend$mean_phi, sp_trend$mean_phi_lcl, sp_trend$mean_phi_ucl),
    sprintf("  Extinction:   %.3f", sp_trend$mean_ext),
    "",
    "DETECTION",
    sprintf("  Mean p = %.3f (%.3f, %.3f)",
            sp_trend$mean_p, sp_trend$mean_p_lcl, sp_trend$mean_p_ucl),
    "",
    "SIGNIFICANT COVARIATE EFFECTS (95% CI excludes 0)"
  )
  
  if (nrow(sig_covs) == 0) {
    txt <- c(txt, "  None")
  } else {
    for (i in 1:nrow(sig_covs)) {
      direction <- ifelse(sig_covs$estimate[i] > 0, "+", "-")
      txt <- c(txt, sprintf("  %s: %.3f (%.3f, %.3f) [%s]",
                             sig_covs$parameter[i], sig_covs$estimate[i],
                             sig_covs$lcl[i], sig_covs$ucl[i], direction))
    }
  }
  
  txt <- c(txt, "",
           "CONVERGENCE",
           sprintf("  Max Rhat: %.3f | Min n.eff: %.0f",
                   sp_trend$max_rhat, sp_trend$min_neff))
  
  writeLines(txt, file.path(FIG_DIR, paste0(sp, "_results.txt")))
  cat("  ", sp, ": figures + text saved\n")
}

# ============================================================================
# MULTI-SPECIES FIGURES (when >1 species loaded)
# ============================================================================

if (nrow(trend_df) > 1) {
  
  cat("\nCreating multi-species figures...\n")
  
  # Detect year range from data
  all_years <- sort(unique(occ_df$year))
  yr_label <- paste0(min(all_years), "\u2013", max(all_years))
  
  # --- All-species occupancy panel ---
  p_panel <- ggplot(occ_df, aes(x = year, y = psi)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2, fill = "#1B9E77") +
    geom_line(color = "#1B9E77", linewidth = 0.8) +
    geom_point(color = "#1B9E77", size = 1.5) +
    facet_wrap(~ species, scales = "free_y", ncol = 4) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = all_years) +
    labs(title = paste0("Occupancy Trends \u2014 All Species (BC, ", yr_label, ")"),
         x = "Year", y = "Occupancy (\u03c8)") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 14),
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  
  n_sp <- length(unique(occ_df$species))
  panel_height <- ceiling(n_sp / 4) * 3 + 1
  ggsave(file.path(FIG_DIR, "ALL_occupancy_panel.png"),
         p_panel, width = 14, height = panel_height, dpi = 300)
  
  # --- Lambda comparison ---
  trend_sorted <- trend_df %>% arrange(lambda_mean)
  trend_sorted$species <- factor(trend_sorted$species, levels = trend_sorted$species)
  trend_sorted$trend_dir <- case_when(
    trend_sorted$lambda_ucl < 1 ~ "Decline",
    trend_sorted$lambda_lcl > 1 ~ "Increase",
    TRUE ~ "Stable"
  )
  
  p_lambda <- ggplot(trend_sorted, aes(x = lambda_mean, y = species, color = trend_dir)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    geom_errorbarh(aes(xmin = lambda_lcl, xmax = lambda_ucl), height = 0.3) +
    geom_point(size = 3) +
    scale_color_manual(values = c("Decline" = "#D73027", "Stable" = "#4575B4",
                                  "Increase" = "#1A9850"), name = "Trend") +
    labs(title = paste0("Occupancy Trends \u2014 All Species (BC, ", yr_label, ")"),
         subtitle = "\u03bb = \u03c8_final / \u03c8_initial | Dashed = stable (\u03bb = 1)",
         x = "Trend (\u03bb)", y = NULL) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(FIG_DIR, "ALL_lambda_comparison.png"),
         p_lambda, width = 10, height = max(4, n_sp * 0.5), dpi = 300)
  
  # --- Covariate effects comparison ---
  coef_df$parameter <- factor(coef_df$parameter, levels = rev(param_order))
  
  p_coef_all <- ggplot(coef_df, aes(x = estimate, y = species, color = significant)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = lcl, xmax = ucl), height = 0.25) +
    geom_point(aes(shape = significant), size = 2) +
    facet_wrap(~ parameter, scales = "free_x", ncol = 3) +
    scale_color_manual(values = c("Yes" = "#D73027", "No" = "gray50"), guide = "none") +
    scale_shape_manual(values = c("Yes" = 16, "No" = 1), guide = "none") +
    labs(title = paste0("Covariate Effects \u2014 All Species (BC, ", yr_label, ")"),
         subtitle = "Filled/red = 95% CI excludes 0",
         x = "Effect (logit scale)", y = NULL) +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 13),
          strip.text = element_text(face = "bold", size = 9),
          axis.text.y = element_text(size = 7))
  
  ggsave(file.path(FIG_DIR, "ALL_covariate_effects.png"),
         p_coef_all, width = 14, height = max(6, n_sp * 0.6), dpi = 300)
  
  # --- Dynamics panel ---
  dyn_years_all <- sort(unique(dyn_df$year_to))
  
  dyn_long <- dyn_df %>%
    select(species, year_to, gamma_mean, gamma_lcl, gamma_ucl,
           phi_mean, phi_lcl, phi_ucl) %>%
    pivot_longer(
      cols = c(gamma_mean, gamma_lcl, gamma_ucl, phi_mean, phi_lcl, phi_ucl),
      names_to = c("rate", ".value"),
      names_pattern = "(gamma|phi)_(mean|lcl|ucl)"
    ) %>%
    mutate(rate = ifelse(rate == "gamma", "Colonization", "Persistence"))
  
  p_dyn_panel <- ggplot(dyn_long, aes(x = year_to, y = mean, color = rate, fill = rate)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.5) +
    facet_wrap(~ species, ncol = 4) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = dyn_years_all) +
    scale_color_manual(values = c("Colonization" = "#D95F02", "Persistence" = "#7570B3")) +
    scale_fill_manual(values = c("Colonization" = "#D95F02", "Persistence" = "#7570B3")) +
    labs(title = paste0("Year-Specific Dynamics \u2014 All Species (BC, ",
                        min(dyn_years_all), "\u2013", max(dyn_years_all), ")"),
         subtitle = "At average covariate values",
         x = "Year", y = "Probability", color = NULL, fill = NULL) +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 13),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  
  ggsave(file.path(FIG_DIR, "ALL_dynamics_panel.png"),
         p_dyn_panel, width = 14, height = panel_height, dpi = 300)
  
  # --- Summary table ---
  summary_table <- trend_df %>%
    mutate(
      lambda_str = sprintf("%.2f (%.2f, %.2f)", lambda_mean, lambda_lcl, lambda_ucl),
      trend = case_when(
        lambda_ucl < 1 ~ "Decline",
        lambda_lcl > 1 ~ "Increase",
        TRUE ~ "Stable"
      )
    ) %>%
    select(species, year_range, lambda_str, trend,
           mean_gamma, mean_phi, mean_p, max_rhat) %>%
    arrange(species)
  
  write.csv(summary_table, file.path(FIG_DIR, "ALL_species_summary.csv"), row.names = FALSE)
  cat("Summary table saved\n")
}

# ============================================================================
# DONE
# ============================================================================


