################################################################################
# BC Bat Monitoring - Visualization for Power Analysis Results (v5)
# 
# Updated to work with 2017-2024 data and v5 covariate structure
#
# Covariates:
#   ψ: Water Nearby (binary) + Dist to Road
#   γ/φ: Dist to Road + Dist to Harvest
#
# RUN: Change SPECIES and run entire script
################################################################################

library(tidyverse)
library(cowplot)
library(scales)

# ============================================================================
# SETTINGS - CHANGE SPECIES HERE
# ============================================================================

SPECIES <- "MYSE"  # Options: LACI, LANO, LABO, MYLU, MYSE

INPUT_DIR <- paste0("outputs/bc_only_2017start/power_analysis/", SPECIES)
OUTPUT_DIR <- INPUT_DIR

# Scenario ordering: Best (Full) to Most Reduced (5-Year)
scenarios <- c("full", "2year", "3year", "4year", "5year")
scenario_labels <- c("Full", "2-Year", "3-Year", "4-Year", "5-Year")
names(scenario_labels) <- scenarios

# Years in study (2017-2024, 8 years)
years <- 2017:2024
nyear <- length(years)

# Years for dynamics (colonization/persistence estimated from year 2 onward)
dynamics_years <- 2018:2024

# ============================================================================
# COLOR PALETTE
# ============================================================================

scenario_colors <- c(
  "full"   = "#2C7BB6",  # Deep blue - reference/best
  "2year"  = "#ABD9E9",  # Light blue
  "3year"  = "#F0E442",  # Yellow
  "4year"  = "#FDAE61",  # Orange
  "5year"  = "#D7191C"   # Red - most reduced
)

COLORS <- scenario_colors

# ============================================================================
# LOAD MODEL FITS
# ============================================================================

cat("Loading model fits for", SPECIES, "...\n")

results_list <- list()
occ_list <- list()
dyn_list <- list()
coef_list <- list()

# Load results CSV if available
results_csv <- paste0(INPUT_DIR, "/results.csv")
if (file.exists(results_csv)) {
  results_from_csv <- read_csv(results_csv, show_col_types = FALSE)
  cat("  Loaded results.csv\n")
} else {
  results_from_csv <- NULL
  cat("  WARNING: results.csv not found\n")
}

for (scenario in scenarios) {
  
  fit_file <- paste0(INPUT_DIR, "/fit_", scenario, ".rds")
  
  if (!file.exists(fit_file)) {
    cat("  WARNING: Missing", fit_file, "\n")
    next
  }
  
  fit <- readRDS(fit_file)
  summ <- fit$summary
  
  cat("  Loaded:", scenario_labels[scenario], "\n")
  
  # Get info from CSV if available
  if (!is.null(results_from_csv)) {
    sc_info <- results_from_csv %>% filter(scenario == !!scenario_labels[scenario])
    if (nrow(sc_info) > 0) {
      n_site_years <- sc_info$site_years[1]
      n_detections <- sc_info$detections[1]
    } else {
      n_site_years <- NA
      n_detections <- NA
    }
  } else {
    n_site_years <- NA
    n_detections <- NA
  }
  
  # Store summary results
  results_list[[scenario]] <- tibble(
    scenario = scenario,
    scenario_label = scenario_labels[scenario],
    n_site_years = n_site_years,
    n_detections = n_detections,
    lambda_mean = summ["lambda", "mean"],
    lambda_lcl = summ["lambda", "2.5%"],
    lambda_ucl = summ["lambda", "97.5%"],
    mean_gamma = summ["mean.gamma", "mean"],
    mean_phi = summ["mean.phi", "mean"],
    mean_p = summ["mean.p", "mean"],
    Rhat_max = max(summ[, "Rhat"], na.rm = TRUE)
  )
  
  # Store annual occupancy
  occ_list[[scenario]] <- tibble(
    scenario = scenario,
    scenario_label = scenario_labels[scenario],
    year = years,
    psi_mean = summ[paste0("psi.fs[", 1:nyear, "]"), "mean"],
    psi_lcl = summ[paste0("psi.fs[", 1:nyear, "]"), "2.5%"],
    psi_ucl = summ[paste0("psi.fs[", 1:nyear, "]"), "97.5%"]
  )
  
  # Store annual dynamics
  dyn_list[[scenario]] <- tibble(
    scenario = scenario,
    scenario_label = scenario_labels[scenario],
    year = dynamics_years,
    gamma_mean = summ[paste0("gamma.avg[", 1:(nyear-1), "]"), "mean"],
    gamma_lcl = summ[paste0("gamma.avg[", 1:(nyear-1), "]"), "2.5%"],
    gamma_ucl = summ[paste0("gamma.avg[", 1:(nyear-1), "]"), "97.5%"],
    phi_mean = summ[paste0("phi.avg[", 1:(nyear-1), "]"), "mean"],
    phi_lcl = summ[paste0("phi.avg[", 1:(nyear-1), "]"), "2.5%"],
    phi_ucl = summ[paste0("phi.avg[", 1:(nyear-1), "]"), "97.5%"]
  )
  
  # Store covariate effects (UPDATED for v5 structure)
  coef_list[[scenario]] <- tibble(
    scenario = scenario,
    scenario_label = scenario_labels[scenario],
    parameter = c(
      "psi_intercept", "psi_water", "psi_road",
      "p_intercept", "p_clutter", "p_temp", "p_julian",
      "col_road", "col_harvest",
      "per_road", "per_harvest"
    ),
    param_label = c(
      "ψ: Intercept", "ψ: Water Nearby", "ψ: Dist to Road",
      "p: Intercept", "p: Clutter", "p: Temperature", "p: Julian",
      "γ: Dist to Road", "γ: Dist to Harvest",
      "φ: Dist to Road", "φ: Dist to Harvest"
    ),
    category = c(
      rep("Initial Occupancy", 3),
      rep("Detection", 4),
      rep("Colonization", 2),
      rep("Persistence", 2)
    ),
    estimate = c(
      summ["beta.psi[1]", "mean"], summ["beta.psi[2]", "mean"], summ["beta.psi[3]", "mean"],
      summ["beta.p[1]", "mean"], summ["beta.p[2]", "mean"], summ["beta.p[3]", "mean"], summ["beta.p[4]", "mean"],
      summ["beta.col.road", "mean"], summ["beta.col.harvest", "mean"],
      summ["beta.per.road", "mean"], summ["beta.per.harvest", "mean"]
    ),
    lcl = c(
      summ["beta.psi[1]", "2.5%"], summ["beta.psi[2]", "2.5%"], summ["beta.psi[3]", "2.5%"],
      summ["beta.p[1]", "2.5%"], summ["beta.p[2]", "2.5%"], summ["beta.p[3]", "2.5%"], summ["beta.p[4]", "2.5%"],
      summ["beta.col.road", "2.5%"], summ["beta.col.harvest", "2.5%"],
      summ["beta.per.road", "2.5%"], summ["beta.per.harvest", "2.5%"]
    ),
    ucl = c(
      summ["beta.psi[1]", "97.5%"], summ["beta.psi[2]", "97.5%"], summ["beta.psi[3]", "97.5%"],
      summ["beta.p[1]", "97.5%"], summ["beta.p[2]", "97.5%"], summ["beta.p[3]", "97.5%"], summ["beta.p[4]", "97.5%"],
      summ["beta.col.road", "97.5%"], summ["beta.col.harvest", "97.5%"],
      summ["beta.per.road", "97.5%"], summ["beta.per.harvest", "97.5%"]
    )
  )
}

# Combine all results
results_df <- bind_rows(results_list)
occ_df <- bind_rows(occ_list)
dyn_df <- bind_rows(dyn_list)
coef_df <- bind_rows(coef_list)

cat("\nLoaded", nrow(results_df), "scenarios\n")

# ============================================================================
# SET FACTOR ORDERING
# ============================================================================

results_df$scenario <- factor(results_df$scenario, levels = scenarios)
results_df$scenario_label <- factor(results_df$scenario_label, levels = scenario_labels)

occ_df$scenario <- factor(occ_df$scenario, levels = scenarios)
occ_df$scenario_label <- factor(occ_df$scenario_label, levels = scenario_labels)

dyn_df$scenario <- factor(dyn_df$scenario, levels = scenarios)
dyn_df$scenario_label <- factor(dyn_df$scenario_label, levels = scenario_labels)

coef_df$scenario <- factor(coef_df$scenario, levels = scenarios)
coef_df$scenario_label <- factor(coef_df$scenario_label, levels = scenario_labels)
coef_df$significant <- ifelse(coef_df$lcl > 0 | coef_df$ucl < 0, "Yes", "No")

# ============================================================================
# PLOT 1: OCCUPANCY TRENDS
# ============================================================================

cat("\nCreating plots...\n")
cat("  1. Occupancy trends\n")

p_occ <- ggplot(occ_df, aes(x = year, y = psi_mean, color = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = psi_lcl, ymax = psi_ucl), alpha = 0.12, color = NA) +
  geom_line(aes(linewidth = scenario)) +
  geom_point(aes(size = scenario)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(breaks = years) +
  scale_color_manual(values = COLORS, labels = scenario_labels, name = "Design") +
  scale_fill_manual(values = COLORS, labels = scenario_labels, name = "Design") +
  scale_linewidth_manual(values = c("full" = 1.5, "2year" = 1, "3year" = 1, 
                                    "4year" = 1, "5year" = 1), guide = "none") +
  scale_size_manual(values = c("full" = 3, "2year" = 2, "3year" = 2, 
                               "4year" = 2, "5year" = 2), guide = "none") +
  labs(
    title = paste(SPECIES, "- Occupancy Trends by Sampling Design"),
    subtitle = "Full design shown with thicker line | Shaded regions = 95% CI",
    x = "Year", y = "Occupancy (ψ)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40", size = 11)
  ) +
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))

ggsave(paste0(OUTPUT_DIR, "/fig_01_occupancy_trends.png"), p_occ, 
       width = 12, height = 7, dpi = 300)

# ============================================================================
# PLOT 2: TREND (λ) COMPARISON
# ============================================================================

cat("  2. Lambda comparison\n")

p_lambda <- ggplot(results_df, aes(x = scenario_label, y = lambda_mean, fill = scenario)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_col(alpha = 0.8, width = 0.7) +
  geom_errorbar(aes(ymin = lambda_lcl, ymax = lambda_ucl), 
                width = 0.25, linewidth = 0.8, color = "gray20") +
  geom_text(aes(label = sprintf("%.2f", lambda_mean)), 
            vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = COLORS, guide = "none") +
  scale_y_continuous(limits = c(0, max(results_df$lambda_ucl, na.rm = TRUE) * 1.15),
                     breaks = seq(0, 2, 0.2)) +
  labs(
    title = paste(SPECIES, "- Trend (λ) by Sampling Design"),
    subtitle = expression(lambda == psi[2024] / psi[2017] * " | Dashed line = no change"),
    x = NULL, y = expression(lambda~"(Trend)")
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40", size = 11),
    axis.text.x = element_text(size = 12, face = "bold")
  )

ggsave(paste0(OUTPUT_DIR, "/fig_02_lambda_comparison.png"), p_lambda, 
       width = 10, height = 6, dpi = 300)

# ============================================================================
# PLOT 3: PRECISION COMPARISON
# ============================================================================

cat("  3. Precision comparison\n")

precision_summary <- results_df %>%
  mutate(
    ci_width = lambda_ucl - lambda_lcl,
    full_ci = ci_width[scenario == "full"],
    precision_loss_pct = round((ci_width / full_ci - 1) * 100, 1)
  )

p_precision <- ggplot(precision_summary, aes(x = scenario_label, y = ci_width, fill = scenario)) +
  geom_col(alpha = 0.85, width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", ci_width)), vjust = -0.5, size = 4) +
  geom_text(aes(label = ifelse(scenario == "full", "", sprintf("+%.0f%%", precision_loss_pct))),
            vjust = 1.5, size = 3.5, color = "black", fontface = "bold") +
  scale_fill_manual(values = COLORS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = paste(SPECIES, "- Trend Estimate Precision by Sampling Design"),
    subtitle = "95% CI width for λ | Smaller = more precise | % shows loss relative to Full",
    x = NULL, y = "95% CI Width"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40", size = 11),
    axis.text.x = element_text(size = 12, face = "bold")
  )

ggsave(paste0(OUTPUT_DIR, "/fig_03_precision_comparison.png"), p_precision, 
       width = 10, height = 6, dpi = 300)

# ============================================================================
# PLOT 4: DYNAMICS - Colonization & Persistence
# ============================================================================

cat("  4. Population dynamics\n")

dyn_long <- dyn_df %>%
  pivot_longer(
    cols = c(gamma_mean, gamma_lcl, gamma_ucl, phi_mean, phi_lcl, phi_ucl),
    names_to = c("rate", "stat"),
    names_pattern = "(gamma|phi)_(mean|lcl|ucl)"
  ) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(
    rate_label = factor(
      ifelse(rate == "gamma", "Colonization (γ)", "Persistence (φ)"),
      levels = c("Colonization (γ)", "Persistence (φ)")
    )
  )

p_dyn <- ggplot(dyn_long, aes(x = year, y = mean, color = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.12, color = NA) +
  geom_line(aes(linewidth = scenario)) +
  geom_point(aes(size = scenario)) +
  facet_wrap(~rate_label, ncol = 1, scales = "free_y") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(breaks = dynamics_years) +
  scale_color_manual(values = COLORS, labels = scenario_labels, name = "Design") +
  scale_fill_manual(values = COLORS, labels = scenario_labels, name = "Design") +
  scale_linewidth_manual(values = c("full" = 1.4, "2year" = 0.9, "3year" = 0.9, 
                                    "4year" = 0.9, "5year" = 0.9), guide = "none") +
  scale_size_manual(values = c("full" = 2.5, "2year" = 1.8, "3year" = 1.8, 
                               "4year" = 1.8, "5year" = 1.8), guide = "none") +
  labs(
    title = paste(SPECIES, "- Population Dynamics by Sampling Design"),
    subtitle = "At average covariate values | Years show transition endpoints",
    x = "Year", y = "Probability"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40", size = 11),
    strip.background = element_rect(fill = "gray92", color = NA),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))

ggsave(paste0(OUTPUT_DIR, "/fig_04_dynamics.png"), p_dyn, 
       width = 12, height = 9, dpi = 300)

# ============================================================================
# PLOT 5: DISTANCE COVARIATE EFFECTS
# ============================================================================

cat("  5. Distance effects\n")

dist_coef <- coef_df %>%
  filter(category %in% c("Colonization", "Persistence"))

dist_coef$param_label <- factor(
  dist_coef$param_label,
  levels = c("φ: Dist to Harvest", "φ: Dist to Road", 
             "γ: Dist to Harvest", "γ: Dist to Road")
)

p_dist <- ggplot(dist_coef, aes(x = estimate, y = param_label, color = scenario)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = lcl, xmax = ucl), height = 0.35, 
                 position = position_dodge(width = 0.7), linewidth = 0.9) +
  geom_point(aes(shape = significant), size = 3.5, 
             position = position_dodge(width = 0.7)) +
  scale_color_manual(values = COLORS, labels = scenario_labels, name = "Design") +
  scale_shape_manual(values = c("Yes" = 16, "No" = 1), 
                     name = "Significant", labels = c("No", "Yes")) +
  labs(
    title = paste(SPECIES, "- Distance Covariate Effects by Sampling Design"),
    subtitle = "Effects on colonization (γ) and persistence (φ) | Filled = CI excludes 0",
    x = "Effect (logit scale)", y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40", size = 11)
  ) +
  guides(color = guide_legend(nrow = 1, order = 1), 
         shape = guide_legend(nrow = 1, order = 2))

ggsave(paste0(OUTPUT_DIR, "/fig_05_distance_effects.png"), p_dist, 
       width = 11, height = 6, dpi = 300)

# ============================================================================
# PLOT 6: ALL COVARIATE EFFECTS
# ============================================================================

cat("  6. All covariate effects\n")

param_order <- c(
  "ψ: Intercept", "ψ: Water Nearby", "ψ: Dist to Road",
  "γ: Dist to Road", "γ: Dist to Harvest",
  "φ: Dist to Road", "φ: Dist to Harvest",
  "p: Intercept", "p: Clutter", "p: Temperature", "p: Julian"
)
coef_df$param_label <- factor(coef_df$param_label, levels = rev(param_order))

p_coef <- ggplot(coef_df, aes(x = estimate, y = param_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = lcl, xmax = ucl, color = scenario), height = 0.4, 
                 position = position_dodge(width = 0.7), linewidth = 0.7) +
  geom_point(aes(color = scenario, shape = significant), size = 2.5, 
             position = position_dodge(width = 0.7)) +
  scale_color_manual(values = COLORS, labels = scenario_labels, name = "Design") +
  scale_shape_manual(values = c("Yes" = 16, "No" = 1), 
                     name = "Significant", labels = c("No", "Yes")) +
  labs(
    title = paste(SPECIES, "- All Covariate Effects by Sampling Design"),
    subtitle = "95% credible intervals | Filled points = CI excludes 0",
    x = "Effect (logit scale)", y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "gray40", size = 10)
  ) +
  guides(color = guide_legend(nrow = 1, order = 1), 
         shape = guide_legend(nrow = 1, order = 2))

ggsave(paste0(OUTPUT_DIR, "/fig_06_all_covariates.png"), p_coef, 
       width = 12, height = 9, dpi = 300)

# ============================================================================
# PLOT 7: PARAMETER PRECISION
# ============================================================================

cat("  7. Parameter precision\n")

precision_params <- coef_df %>%
  mutate(ci_width = ucl - lcl) %>%
  filter(parameter %in% c("psi_intercept", "col_road", "col_harvest", 
                          "per_road", "per_harvest", "p_intercept")) %>%
  mutate(
    param_short = case_when(
      parameter == "psi_intercept" ~ "Initial Occ. (ψ₀)",
      parameter == "col_road" ~ "γ: Road",
      parameter == "col_harvest" ~ "γ: Harvest",
      parameter == "per_road" ~ "φ: Road",
      parameter == "per_harvest" ~ "φ: Harvest",
      parameter == "p_intercept" ~ "Detection (p₀)"
    )
  )

precision_params$param_short <- factor(
  precision_params$param_short,
  levels = c("Initial Occ. (ψ₀)", "γ: Road", "γ: Harvest", 
             "φ: Road", "φ: Harvest", "Detection (p₀)")
)

p_param_prec <- ggplot(precision_params, 
                       aes(x = scenario_label, y = ci_width, fill = scenario)) +
  geom_col(width = 0.75, alpha = 0.85) +
  facet_wrap(~param_short, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = COLORS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = paste(SPECIES, "- Parameter Precision by Sampling Design"),
    subtitle = "95% CI width | Smaller = more precise",
    x = NULL, y = "CI Width"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "gray40", size = 10),
    strip.background = element_rect(fill = "gray92", color = NA),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
  )

ggsave(paste0(OUTPUT_DIR, "/fig_07_parameter_precision.png"), p_param_prec, 
       width = 12, height = 7, dpi = 300)

# ============================================================================
# PLOT 8: SUMMARY DASHBOARD
# ============================================================================

cat("  8. Summary dashboard\n")

p_dash_lambda <- ggplot(results_df, aes(x = scenario_label, y = lambda_mean, fill = scenario)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_col(alpha = 0.8, width = 0.7) +
  geom_errorbar(aes(ymin = lambda_lcl, ymax = lambda_ucl), width = 0.2, color = "gray30") +
  scale_fill_manual(values = COLORS, guide = "none") +
  labs(title = "A. Trend (λ)", x = NULL, y = "λ") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

p_dash_prec <- ggplot(precision_summary, aes(x = scenario_label, y = ci_width, fill = scenario)) +
  geom_col(alpha = 0.8, width = 0.7) +
  scale_fill_manual(values = COLORS, guide = "none") +
  labs(title = "B. Precision (λ CI Width)", x = NULL, y = "CI Width") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

sample_size_df <- results_df %>%
  mutate(pct_full = n_site_years / max(n_site_years, na.rm = TRUE) * 100)

p_dash_sample <- ggplot(sample_size_df, aes(x = scenario_label, y = n_site_years, fill = scenario)) +
  geom_col(alpha = 0.8, width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", pct_full)), vjust = -0.3, size = 3) +
  scale_fill_manual(values = COLORS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "C. Site-Years Sampled", x = NULL, y = "N") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

p_dashboard <- plot_grid(p_dash_lambda, p_dash_prec, p_dash_sample, nrow = 1)

title <- ggdraw() + 
  draw_label(paste(SPECIES, "- Sampling Design Comparison Summary (2017-2024)"), 
             fontface = "bold", size = 14)

p_dashboard_final <- plot_grid(title, p_dashboard, ncol = 1, rel_heights = c(0.1, 1))

ggsave(paste0(OUTPUT_DIR, "/fig_08_summary_dashboard.png"), p_dashboard_final, 
       width = 12, height = 5, dpi = 300)

# ============================================================================
# FULL MODEL SUMMARY PLOTS
# ============================================================================

cat("  9-11. Full model summary plots\n")

full_fit <- readRDS(paste0(INPUT_DIR, "/fit_full.rds"))
full_summ <- full_fit$summary

# Full model occupancy
psi_full_df <- data.frame(
  year = years,
  psi = full_summ[paste0("psi.fs[", 1:nyear, "]"), "mean"],
  lcl = full_summ[paste0("psi.fs[", 1:nyear, "]"), "2.5%"],
  ucl = full_summ[paste0("psi.fs[", 1:nyear, "]"), "97.5%"]
)

p_full_occ <- ggplot(psi_full_df, aes(x = year, y = psi)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.25, fill = "#1B9E77") +
  geom_line(linewidth = 1.3, color = "#1B9E77") +
  geom_point(size = 3.5, color = "#1B9E77") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_x_continuous(breaks = years) +
  labs(
    title = paste(SPECIES, "- Occupancy Trend (Full Design)"),
    subtitle = sprintf("λ = %.2f (%.2f, %.2f)", 
                       full_summ["lambda", "mean"],
                       full_summ["lambda", "2.5%"],
                       full_summ["lambda", "97.5%"]),
    x = "Year", y = "Occupancy (ψ)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12, color = "gray30")
  )

ggsave(paste0(OUTPUT_DIR, "/fig_full_01_occupancy.png"), p_full_occ, 
       width = 10, height = 6, dpi = 300)

# Full model dynamics
dyn_full_df <- data.frame(
  year = dynamics_years,
  col_mean = full_summ[paste0("gamma.avg[", 1:(nyear-1), "]"), "mean"],
  col_lcl = full_summ[paste0("gamma.avg[", 1:(nyear-1), "]"), "2.5%"],
  col_ucl = full_summ[paste0("gamma.avg[", 1:(nyear-1), "]"), "97.5%"],
  per_mean = full_summ[paste0("phi.avg[", 1:(nyear-1), "]"), "mean"],
  per_lcl = full_summ[paste0("phi.avg[", 1:(nyear-1), "]"), "2.5%"],
  per_ucl = full_summ[paste0("phi.avg[", 1:(nyear-1), "]"), "97.5%"]
)

p_full_dyn <- ggplot(dyn_full_df) +
  geom_ribbon(aes(x = year, ymin = col_lcl, ymax = col_ucl), alpha = 0.2, fill = "#D95F02") +
  geom_line(aes(x = year, y = col_mean, color = "Colonization"), linewidth = 1.2) +
  geom_point(aes(x = year, y = col_mean, color = "Colonization"), size = 3) +
  geom_ribbon(aes(x = year, ymin = per_lcl, ymax = per_ucl), alpha = 0.2, fill = "#7570B3") +
  geom_line(aes(x = year, y = per_mean, color = "Persistence"), linewidth = 1.2) +
  geom_point(aes(x = year, y = per_mean, color = "Persistence"), size = 3) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_x_continuous(breaks = dynamics_years) +
  scale_color_manual(values = c("Colonization" = "#D95F02", "Persistence" = "#7570B3")) +
  labs(
    title = paste(SPECIES, "- Population Dynamics (Full Design)"),
    subtitle = "At average covariate values | Years show transition endpoints",
    x = "Year", y = "Probability", color = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30")
  )

ggsave(paste0(OUTPUT_DIR, "/fig_full_02_dynamics.png"), p_full_dyn, 
       width = 10, height = 6, dpi = 300)

# Full model covariates
full_coef_df <- data.frame(
  parameter = c(
    "ψ: Intercept", "ψ: Water Nearby", "ψ: Dist to Road",
    "p: Intercept", "p: Clutter", "p: Temperature", "p: Julian",
    "γ: Dist to Road", "γ: Dist to Harvest",
    "φ: Dist to Road", "φ: Dist to Harvest"
  ),
  estimate = c(
    full_summ["beta.psi[1]", "mean"], full_summ["beta.psi[2]", "mean"], full_summ["beta.psi[3]", "mean"],
    full_summ["beta.p[1]", "mean"], full_summ["beta.p[2]", "mean"], full_summ["beta.p[3]", "mean"], full_summ["beta.p[4]", "mean"],
    full_summ["beta.col.road", "mean"], full_summ["beta.col.harvest", "mean"],
    full_summ["beta.per.road", "mean"], full_summ["beta.per.harvest", "mean"]
  ),
  lcl = c(
    full_summ["beta.psi[1]", "2.5%"], full_summ["beta.psi[2]", "2.5%"], full_summ["beta.psi[3]", "2.5%"],
    full_summ["beta.p[1]", "2.5%"], full_summ["beta.p[2]", "2.5%"], full_summ["beta.p[3]", "2.5%"], full_summ["beta.p[4]", "2.5%"],
    full_summ["beta.col.road", "2.5%"], full_summ["beta.col.harvest", "2.5%"],
    full_summ["beta.per.road", "2.5%"], full_summ["beta.per.harvest", "2.5%"]
  ),
  ucl = c(
    full_summ["beta.psi[1]", "97.5%"], full_summ["beta.psi[2]", "97.5%"], full_summ["beta.psi[3]", "97.5%"],
    full_summ["beta.p[1]", "97.5%"], full_summ["beta.p[2]", "97.5%"], full_summ["beta.p[3]", "97.5%"], full_summ["beta.p[4]", "97.5%"],
    full_summ["beta.col.road", "97.5%"], full_summ["beta.col.harvest", "97.5%"],
    full_summ["beta.per.road", "97.5%"], full_summ["beta.per.harvest", "97.5%"]
  ),
  category = c(
    rep("Initial Occupancy", 3), rep("Detection", 4),
    rep("Colonization", 2), rep("Persistence", 2)
  )
)

full_coef_df$parameter <- factor(full_coef_df$parameter, levels = rev(param_order))
full_coef_df$significant <- ifelse(full_coef_df$lcl > 0 | full_coef_df$ucl < 0, "Yes", "No")

category_colors <- c(
  "Initial Occupancy" = "#1B9E77", "Detection" = "#666666",
  "Colonization" = "#D95F02", "Persistence" = "#7570B3"
)

p_full_coef <- ggplot(full_coef_df, aes(x = estimate, y = parameter, color = category)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = lcl, xmax = ucl), height = 0.35, linewidth = 0.9) +
  geom_point(aes(shape = significant), size = 4) +
  scale_color_manual(values = category_colors, name = "Process") +
  scale_shape_manual(values = c("Yes" = 16, "No" = 1), name = "Significant") +
  labs(
    title = paste(SPECIES, "- Covariate Effects (Full Design)"),
    subtitle = "95% credible intervals | Filled = CI excludes 0",
    x = "Effect (logit scale)", y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30")
  )

ggsave(paste0(OUTPUT_DIR, "/fig_full_03_covariates.png"), p_full_coef, 
       width = 11, height = 8, dpi = 300)

# ============================================================================
# PRINT SUMMARY
# ============================================================================

cat("\n", rep("=", 68), "\n", sep = "")
cat("RESULTS SUMMARY:", SPECIES, "(2017-2024)\n")
cat(rep("=", 68), "\n\n", sep = "")

summary_print <- results_df %>%
  mutate(
    ci_width = round(lambda_ucl - lambda_lcl, 3),
    lambda_str = sprintf("%.2f (%.2f, %.2f)", lambda_mean, lambda_lcl, lambda_ucl)
  ) %>%
  select(scenario_label, n_site_years, lambda_str, ci_width)

print(summary_print)

cat("\nPrecision loss (CI width relative to Full):\n")
print(precision_summary %>% select(scenario_label, ci_width, precision_loss_pct))

cat("\n", rep("=", 68), "\n", sep = "")
cat("COMPLETE! Figures saved to:", OUTPUT_DIR, "\n")
cat(rep("=", 68), "\n", sep = "")

