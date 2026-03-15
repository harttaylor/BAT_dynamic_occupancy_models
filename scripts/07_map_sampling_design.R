################################################################################
# BC Bat Monitoring - Rotation Schedule Maps (IMPROVED)
# 
# Uses BCRegions.shp for accurate region boundaries
# IMPROVEMENTS:
# - Muted/pastel region colors (not distracting)
# - Better terminology: "Rotation Schedule" not "Panel"
# - Clearer site symbols
# - Core vs Rotating terminology
#
################################################################################

library(tidyverse)
library(sf)
library(cowplot)

# ============================================================================
# CONFIGURATION
# ============================================================================

VALID_QUADRANTS <- c("NE", "NW", "SE", "SW")
OUTPUT_DIR <- "outputs/bc_only/maps"

# Site classification from ownership CSV
ROTATING_GRTS <- c(
  186490, 137274, 77162, 243562, 226090, 128058, 114602, 143274,
  16426, 175354, 37274, 283586, 212906, 6202, 172090, 986, 
  265749, 317428, 38442, 132650, 115322, 206394, 26170, 
  230186, 17706, 145578, 58842, 110378
)

ALWAYS_ON_GRTS <- c(
  47210, 89962, 20330, 271731, 241002, 151210, 216746, 273050, 
  195242, 19514, 219626, 252394, 127850, 47322, 294832, 269777, 
  20730, 9466, 42026, 53210, 186330, 1194, 79018, 218282, 5994, 
  73002, 147306, 12522, 86074, 326072, 57322, 73786, 60586, 
  70250, 268521, 93290, 143338, 212266, 194666, 159722, 111466
)

# ============================================================================
# MUTED COLOR PALETTE (not distracting)
# ============================================================================

# Soft, muted region colors - visible but not overwhelming
region_colors_muted <- c(
  "Northern" = "#C7E9C0",      # Soft sage green
  
  "Southcentral" = "#DADAEB",  # Soft lavender
  "Southeast" = "#FDD0A2",    # Soft peach
  "Southwest" = "#C6DBEF"     # Soft sky blue
)

# Site colors - high contrast for visibility
site_colors <- list(
  core = "#2166AC",           # Strong blue - always monitored
  rotating_on = "#D95F02",    # Orange - rotating & sampled
  rotating_off = "white"      # White with gray border - not sampled
)

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading data...\n")

# Load BC regions shapefile
bc_regions <- st_read("data/raw/BCRegions.shp", quiet = TRUE)
cat("Loaded BC regions shapefile with", nrow(bc_regions), "regions\n")

# Load site data
raw_data <- read.csv("data/BatActivity_2016to2024.csv")

site_info <- raw_data %>%
  filter(Region != "Alaska", !is.na(Region)) %>%
  mutate(
    site = paste0(GRTS.Cell.ID, "_", Quadrant),
    grts_cell = GRTS.Cell.ID
  ) %>%
  filter(Quadrant %in% VALID_QUADRANTS) %>%
  group_by(site) %>%
  summarize(
    grts_cell = first(grts_cell),
    region = first(na.omit(Region)),
    lat = mean(Lat, na.rm = TRUE),
    long = mean(Long, na.rm = TRUE),
    n_years = n_distinct(as.integer(format(as.Date(Night), "%Y"))),
    .groups = "drop"
  ) %>%
  filter(!is.na(lat), !is.na(long), n_years >= 2) %>%
  mutate(
    site_type = case_when(
      grts_cell %in% ROTATING_GRTS ~ "Rotating",
      grts_cell %in% ALWAYS_ON_GRTS ~ "Core",
      TRUE ~ "Core"  # Default to core
    )
  )

cat("Loaded", nrow(site_info), "BC sites\n")
cat("  Core (always monitored):", sum(site_info$site_type == "Core"), "\n")
cat("  Rotating:", sum(site_info$site_type == "Rotating"), "\n")

# Assign rotation groups within region
set.seed(42)
site_info <- site_info %>%
  group_by(region, site_type) %>%
  mutate(
    random_order = sample(n()),
    rotation_group = ifelse(site_type == "Core", 0, ((random_order - 1) %% 5) + 1)
  ) %>%
  ungroup() %>%
  dplyr::select(-random_order)

# Convert to sf
sites_sf <- site_info %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326)

# ============================================================================
# MAP FUNCTIONS
# ============================================================================

# Terminology mapping
schedule_names <- c(
  "1" = "Annual",
  "2" = "2-Year Rotation",
  "3" = "3-Year Rotation",
  "4" = "4-Year Rotation",
  "5" = "5-Year Rotation"
)

# Get active rotation groups for a schedule/year
get_active_groups <- function(year_in_cycle, n_year_rotation) {
  # n_year_rotation: 1 = annual, 2 = 2-year, etc.
  # Returns which rotation groups (1-5) are sampled
  
  if (n_year_rotation == 1) return(1:5)  # Annual: all rotating sites
  
  if (n_year_rotation == 2) {
    if (year_in_cycle == 1) return(c(1, 2, 3))
    else return(c(4, 5))
  } else if (n_year_rotation == 3) {
    if (year_in_cycle == 1) return(c(1, 2))
    else if (year_in_cycle == 2) return(c(3, 4))
    else return(c(5))
  } else if (n_year_rotation == 4) {
    if (year_in_cycle == 1) return(c(1, 2))
    else if (year_in_cycle == 2) return(c(3))
    else if (year_in_cycle == 3) return(c(4))
    else return(c(5))
  } else if (n_year_rotation == 5) {
    return(year_in_cycle)  # Only one group per year
  }
  return(1:5)
}

make_schedule_map <- function(n_year_rotation, year_in_cycle, title = NULL, 
                              show_legend = FALSE, compact = FALSE) {
  
  # Determine which groups are active
  if (n_year_rotation == 1) {
    # Annual: all sites sampled
    active_groups <- 1:5
  } else {
    active_groups <- get_active_groups(year_in_cycle, n_year_rotation)
  }
  
  # Prepare site data
  sites_plot <- sites_sf %>%
    mutate(
      is_sampled = (site_type == "Core") | (rotation_group %in% active_groups),
      status = case_when(
        site_type == "Core" ~ "Core (always monitored)",
        is_sampled ~ "Rotating (sampled this year)",
        TRUE ~ "Rotating (not sampled)"
      ),
      status = factor(status, levels = c("Core (always monitored)", 
                                         "Rotating (sampled this year)",
                                         "Rotating (not sampled)"))
    )
  
  # Count sites
  n_core <- sum(sites_plot$site_type == "Core")
  n_rotating_on <- sum(sites_plot$is_sampled & sites_plot$site_type == "Rotating")
  n_total <- n_core + n_rotating_on
  
  # Generate title if not provided
  if (is.null(title)) {
    if (n_year_rotation == 1) {
      title <- "Annual (All Sites)"
    } else {
      title <- paste0(n_year_rotation, "-Year (Year ", year_in_cycle, ")")
    }
  }
  
  subtitle <- sprintf("%d core + %d rotating = %d sites", n_core, n_rotating_on, n_total)
  
  # Base plot
  p <- ggplot() +
    # Region polygons with muted colors
    geom_sf(data = bc_regions, aes(fill = Region), 
            color = "gray60", linewidth = 0.4, alpha = 0.7) +
    scale_fill_manual(values = region_colors_muted, guide = "none") +
    
    # Sites not sampled (gray border, white fill)
    geom_sf(data = filter(sites_plot, status == "Rotating (not sampled)"),
            fill = "white", color = "gray50", size = 2, shape = 21, stroke = 0.5) +
    
    # Core sites (blue)
    geom_sf(data = filter(sites_plot, status == "Core (always monitored)"),
            fill = site_colors$core, color = "black", size = 2.5, shape = 21, stroke = 0.6) +
    
    # Rotating sites sampled (orange)
    geom_sf(data = filter(sites_plot, status == "Rotating (sampled this year)"),
            fill = site_colors$rotating_on, color = "black", size = 2.5, shape = 21, stroke = 0.6) +
    
    # Zoom to BC
    coord_sf(xlim = c(-136, -114), ylim = c(48.3, 60), expand = FALSE)
  
  # Styling
  if (compact) {
    p <- p +
      labs(title = title, subtitle = subtitle) +
      theme_void(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
        plot.subtitle = element_text(size = 8, hjust = 0.5, color = "gray40"),
        plot.margin = margin(2, 2, 2, 2)
      )
  } else {
    p <- p +
      labs(title = title, subtitle = subtitle) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
        panel.grid = element_line(color = "gray90", linewidth = 0.2),
        panel.background = element_rect(fill = "#F8F8F8", color = NA),
        axis.text = element_text(size = 8, color = "gray50"),
        axis.title = element_blank()
      )
  }
  
  return(p)
}

# ============================================================================
# CREATE FIGURES
# ============================================================================

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- LEGEND ---
legend_text <- "●Core (always monitored)    ●Rotating (sampled)    ○Rotating (not sampled)"

create_legend_row <- function() {
  ggdraw() +
    draw_label(
      "● Core sites (always monitored)    ● Rotating sites (sampled)    ○ Rotating sites (not sampled this year)",
      size = 11, color = "gray30"
    )
}

# Add colored legend
create_colored_legend <- function() {
  legend_df <- data.frame(
    x = c(1, 2, 3),
    status = factor(c("Core (always monitored)", 
                      "Rotating (sampled)", 
                      "Rotating (not sampled)"),
                    levels = c("Core (always monitored)", 
                               "Rotating (sampled)", 
                               "Rotating (not sampled)"))
  )
  
  p <- ggplot(legend_df, aes(x = x, y = 1)) +
    geom_point(aes(fill = status, color = status), shape = 21, size = 5, stroke = 0.8) +
    scale_fill_manual(
      values = c(
        "Core (always monitored)" = site_colors$core,
        "Rotating (sampled)" = site_colors$rotating_on,
        "Rotating (not sampled)" = "white"
      )
    ) +
    scale_color_manual(
      values = c(
        "Core (always monitored)" = "black",
        "Rotating (sampled)" = "black",
        "Rotating (not sampled)" = "gray50"
      )
    ) +
    geom_text(aes(label = status), hjust = -0.2, size = 3.5) +
    scale_x_continuous(limits = c(0.5, 6)) +
    theme_void() +
    theme(legend.position = "none")
  
  return(p)
}

# --- FIGURE 1: SCHEDULE COMPARISON (Year 1 of each) ---
cat("\nCreating schedule comparison figure...\n")

p_annual <- make_schedule_map(1, 1, "Annual\n(All Sites)", compact = TRUE)
p_2yr <- make_schedule_map(2, 1, "2-Year Rotation\n(Year 1 of 2)", compact = TRUE)
p_3yr <- make_schedule_map(3, 1, "3-Year Rotation\n(Year 1 of 3)", compact = TRUE)
p_4yr <- make_schedule_map(4, 1, "4-Year Rotation\n(Year 1 of 4)", compact = TRUE)
p_5yr <- make_schedule_map(5, 1, "5-Year Rotation\n(Year 1 of 5)", compact = TRUE)

comparison_maps <- plot_grid(p_annual, p_2yr, p_3yr, p_4yr, p_5yr,
                             ncol = 5, align = "hv")

title_comparison <- ggdraw() + 
  draw_label("Rotation Schedule Comparison: Spatial Coverage in Year 1", 
             fontface = "bold", size = 16)

legend_comparison <- ggdraw() +
  draw_label("● Core (always monitored)     ● Rotating (sampled this year)     ○ Rotating (not sampled)",
             size = 11, color = "gray30")

fig_comparison <- plot_grid(
  title_comparison,
  comparison_maps,
  legend_comparison,
  ncol = 1,
  rel_heights = c(0.07, 1, 0.05)
)

ggsave(paste0(OUTPUT_DIR, "/01_rotation_schedule_comparison.png"), fig_comparison,
       width = 20, height = 5.5, dpi = 300, bg = "white")
cat("Saved: 01_rotation_schedule_comparison.png\n")

# --- FIGURE 2: 2-YEAR ROTATION DETAIL ---
cat("Creating 2-year rotation detail figure...\n")

p_2yr_y1 <- make_schedule_map(2, 1, "Year 1: Groups 1-3 Sampled", compact = FALSE)
p_2yr_y2 <- make_schedule_map(2, 2, "Year 2: Groups 4-5 Sampled", compact = FALSE)

maps_2yr <- plot_grid(p_2yr_y1, p_2yr_y2, ncol = 2, align = "hv")

title_2yr <- ggdraw() + 
  draw_label("2-Year Rotation Schedule", fontface = "bold", size = 16) 

subtitle_2yr <- ggdraw() +
  draw_label("Sites alternate between sampled and not-sampled years. Core sites monitored every year.",
             size = 11, color = "gray40")

fig_2yr <- plot_grid(
  title_2yr,
  subtitle_2yr,
  maps_2yr,
  legend_comparison,
  ncol = 1,
  rel_heights = c(0.06, 0.04, 1, 0.05)
)

ggsave(paste0(OUTPUT_DIR, "/02_2year_rotation_detail.png"), fig_2yr,
       width = 14, height = 8, dpi = 300, bg = "white")
cat("Saved: 02_2year_rotation_detail.png\n")

# --- FIGURE 3: 5-YEAR ROTATION DETAIL ---
cat("Creating 5-year rotation detail figure...\n")

p5_maps <- list()
for (yr in 1:5) {
  p5_maps[[yr]] <- make_schedule_map(5, yr, paste0("Year ", yr, " of 5"), compact = TRUE)
}

maps_5yr <- plot_grid(plotlist = p5_maps, ncol = 5, align = "hv")

title_5yr <- ggdraw() + 
  draw_label("5-Year Rotation Schedule", fontface = "bold", size = 16)

subtitle_5yr <- ggdraw() +
  draw_label("Each rotating site sampled once every 5 years. Core sites monitored every year.",
             size = 11, color = "gray40")

fig_5yr <- plot_grid(
  title_5yr,
  subtitle_5yr,
  maps_5yr,
  legend_comparison,
  ncol = 1,
  rel_heights = c(0.06, 0.04, 1, 0.05)
)

ggsave(paste0(OUTPUT_DIR, "/03_5year_rotation_detail.png"), fig_5yr,
       width = 20, height = 6, dpi = 300, bg = "white")
cat("Saved: 03_5year_rotation_detail.png\n")

# --- FIGURE 4: ANNUAL SAMPLING (Reference) ---
cat("Creating annual reference figure...\n")

p_annual_full <- ggplot() +
  geom_sf(data = bc_regions, aes(fill = Region), 
          color = "gray60", linewidth = 0.5, alpha = 0.7) +
  scale_fill_manual(values = region_colors_muted, name = "Region") +
  geom_sf(data = filter(sites_sf, site_type == "Core"),
          aes(color = "Core (always monitored)"), 
          fill = site_colors$core, size = 3, shape = 21, stroke = 0.7) +
  geom_sf(data = filter(sites_sf, site_type == "Rotating"),
          aes(color = "Rotating sites"), 
          fill = site_colors$rotating_on, size = 3, shape = 21, stroke = 0.7) +
  scale_color_manual(
    values = c("Core (always monitored)" = "black", "Rotating sites" = "black"),
    name = "Site Type",
    guide = guide_legend(override.aes = list(
      fill = c(site_colors$core, site_colors$rotating_on),
      shape = 21, size = 4
    ))
  ) +
  coord_sf(xlim = c(-136, -114), ylim = c(48.3, 60), expand = FALSE) +
  labs(
    title = "Annual Sampling: All Sites Monitored Every Year",
    subtitle = sprintf("%d core sites + %d rotating sites = %d total",
                       sum(site_info$site_type == "Core"),
                       sum(site_info$site_type == "Rotating"),
                       nrow(site_info))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    legend.position = "right",
    panel.grid = element_line(color = "gray90", linewidth = 0.2),
    panel.background = element_rect(fill = "#F8F8F8", color = NA)
  )

ggsave(paste0(OUTPUT_DIR, "/04_annual_sampling_reference.png"), p_annual_full,
       width = 12, height = 10, dpi = 300, bg = "white")
cat("Saved: 04_annual_sampling_reference.png\n")

# ============================================================================
# SITE SUMMARY
# ============================================================================

cat("\n", rep("=", 60), "\n", sep = "")
cat("SITE CLASSIFICATION SUMMARY\n")
cat(rep("=", 60), "\n\n", sep = "")

cat("Overall:\n")
cat(sprintf("  Core (always monitored): %d sites\n", sum(site_info$site_type == "Core")))
cat(sprintf("  Rotating: %d sites\n", sum(site_info$site_type == "Rotating")))
cat(sprintf("  Total: %d sites\n\n", nrow(site_info)))

cat("By region:\n")
site_info %>%
  count(region, site_type) %>%
  pivot_wider(names_from = site_type, values_from = n, values_fill = 0) %>%
  mutate(Total = Core + Rotating) %>%
  print()

cat("\nRotating sites by group (for 5-year rotation):\n")
site_info %>%
  filter(site_type == "Rotating") %>%
  count(rotation_group) %>%
  print()

# Save site classification
write_csv(site_info, paste0(OUTPUT_DIR, "/site_classification_with_groups.csv"))

cat("\n", rep("=", 60), "\n", sep = "")
cat("MAPS COMPLETE\n")
cat("Output:", OUTPUT_DIR, "\n")
cat(rep("=", 60), "\n", sep = "")