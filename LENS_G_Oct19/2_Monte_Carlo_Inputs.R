#### 10,000 Monte Carlo simulations for each island

# ==============================================================================
# MASTER HYDROLOGICAL VALIDATION ENGINE & PUBLICATION GRAPHICS
# ==============================================================================
library(tidyverse)
library(lubridate)
library(patchwork)
library(scales)

# Set global project workspace root directory
setwd("Q:/Hydro/projects/indopacom_PRIVATE/ILANA") 

# Load the master empirical island data matrix
island_features_raw <- read.csv("Data_Files/final_results.csv")

# Isolate the exact validation target list and recode FSM to Micronesia
island_features2 <- island_features_raw %>%
  filter(Island %in% c(
    "Kwajelein", "Roi-Namur", "Eneu", "Laura", "Bikini", "Enjebi",
    "Falalop", "Deke", "Pingelap", "Kalap", "Ngatik", "Buariki",
    "Buota", "Bonriki", "Betio", "Abaiang", "Butaritari", "Maiana",
    "Marakei", "Nonouti", "Onotoa", "Cantonment", "Home", "South",
    "Hoarafushi", "Ihavandhoo", "Dhidhdhoo", "Kelaa", "Filladhoo", "Baarah",
    "Hanimaadhoo", "Nolhivaranfaru", "Nolhivaram", "Kulhudhuffushi", "Kumundhoo", "Hithadhoo",
    "Maradhoo", "Feydhoo", "Gan (Addu)", "Thinadhoo", "Gan (Laamu)", "Holhudhoo",
    "Velidhoo", "Cooper-Meng", "Eastern", "Wale", "Motu_Kotawa", "Moto_Ko"
  )) %>%
  mutate(Region = recode(Region, "FSM" = "Micronesia"))

# Master Analytical Geometry K-Lookup Grid (Embedded directly to prevent missing file crashes)
island_width_k <- data.frame(
  width = c(150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 
            850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1300, 1400, 1500, 1600,
            1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500),
  k25  = c(1.6, 1.58, 1.54, 1.47, 1.36, 1.24, 1.13, 1.07, 1.03, 1.02, 1.01, 1.01, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
  k50  = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
  k75  = c(0.74, 0.75, 0.77, 0.79, 0.81, 0.83, 0.86, 0.88, 0.91, 0.93, 0.94, 0.96, 0.97, 0.98, 0.98, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
  k100 = c(0.61, 0.62, 0.63, 0.66, 0.68, 0.72, 0.75, 0.79, 0.83, 0.87, 0.9, 0.92, 0.94, 0.95, 0.96, 0.97, 0.98, 0.98, 0.98, 0.98, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99),
  k125 = c(0.52, 0.53, 0.55, 0.57, 0.6, 0.63, 0.67, 0.71, 0.76, 0.8, 0.84, 0.87, 0.9, 0.92, 0.94, 0.95, 0.96, 0.97, 0.97, 0.97, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98),
  k150 = c(0.46, 0.47, 0.49, 0.51, 0.53, 0.57, 0.6, 0.65, 0.69, 0.74, 0.78, 0.82, 0.86, 0.88, 0.91, 0.93, 0.94, 0.95, 0.96, 0.96, 0.97, 0.97, 0.97, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98),
  k175 = c(0.41, 0.42, 0.44, 0.46, 0.48, 0.51, 0.55, 0.59, 0.63, 0.68, 0.72, 0.77, 0.81, 0.84, 0.87, 0.89, 0.91, 0.93, 0.94, 0.95, 0.96, 0.96, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97),
  k200 = c(0.38, 0.39, 0.4, 0.42, 0.44, 0.47, 0.5, 0.53, 0.57, 0.62, 0.66, 0.71, 0.75, 0.79, 0.83, 0.86, 0.88, 0.9, 0.92, 0.93, 0.94, 0.95, 0.96, 0.96, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97),
  k225 = c(0.35, 0.36, 0.37, 0.39, 0.4, 0.43, 0.45, 0.49, 0.53, 0.57, 0.61, 0.66, 0.70, 0.74, 0.78, 0.82, 0.85, 0.87, 0.89, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96),
  k250 = c(0.33, 0.33, 0.34, 0.36, 0.37, 0.4, 0.42, 0.45, 0.48, 0.52, 0.56, 0.61, 0.65, 0.70, 0.74, 0.77, 0.81, 0.84, 0.86, 0.88, 0.90, 0.91, 0.93, 0.94, 0.95, 0.95, 0.95, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96),
  k275 = c(0.31, 0.31, 0.32, 0.34, 0.35, 0.37, 0.39, 0.42, 0.45, 0.48, 0.52, 0.56, 0.61, 0.65, 0.69, 0.73, 0.77, 0.80, 0.83, 0.85, 0.87, 0.89, 0.91, 0.93, 0.94, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.96, 0.96, 0.96),
  k300 = c(0.29, 0.30, 0.31, 0.32, 0.33, 0.35, 0.37, 0.39, 0.42, 0.45, 0.49, 0.52, 0.56, 0.61, 0.65, 0.69, 0.73, 0.76, 0.79, 0.82, 0.85, 0.87, 0.90, 0.92, 0.94, 0.94, 0.95, 0.95, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96),
  k325 = c(0.27, 0.28, 0.29, 0.30, 0.31, 0.33, 0.35, 0.37, 0.40, 0.42, 0.46, 0.49, 0.53, 0.57, 0.61, 0.65, 0.69, 0.73, 0.76, 0.79, 0.82, 0.85, 0.89, 0.92, 0.94, 0.95, 0.96, 0.97, 0.97, 0.97, 0.98, 0.98, 0.98, 0.98, 0.98),
  k400 = c(0.24, 0.25, 0.25, 0.26, 0.27, 0.29, 0.30, 0.32, 0.34, 0.37, 0.39, 0.42, 0.45, 0.49, 0.53, 0.56, 0.60, 0.64, 0.68, 0.71, 0.74, 0.77, 0.82, 0.86, 0.89, 0.91, 0.92, 0.93, 0.94, 0.94, 0.95, 0.95, 0.95, 0.95, 0.95)
)


# ==============================================================================
# 2. CONTINUOUS SCALING MATHEMATICAL FUNCTIONS
# ==============================================================================

# Core Analytical Base Formula
Bailey_Model <- function(Y, Z_TD, R, B, K, C) {
  ((((Z_TD - Y) * R) / (B + R)) + Y) * K * C
}

# Empirical Exponential Fit Equations (Replaces external file lookups)
get_C_factor <- function(width) { 1 - 0.15 * exp(-0.00288 * width) }
get_B_factor <- function(width) { 5 * exp(-0.00461 * width) }
get_Y_factor <- function(width) { -20 * (1 - exp(-0.001386 * width)) }

# 2D Linear Numerical Interpolation Function for K Coefficients
interpolate_K <- function(width, hc) {
  nearest_width <- island_width_k$width[which.min(abs(island_width_k$width - width))]
  row <- island_width_k[island_width_k$width == nearest_width, ]
  k_vals <- as.numeric(row[-1])
  hc_levels <- as.numeric(gsub("k", "", names(row)[-1]))
  
  valid <- which(!is.na(k_vals) & !is.na(hc_levels))
  if (length(valid) < 2) return(NA_real_)
  approx(x = hc_levels[valid], y = k_vals[valid], xout = hc, rule = 2)$y
}


# ==============================================================================
# 3. STOCHASTIC MONTE CARLO CORE ENGINE
# ==============================================================================

MC_Bailey <- function(low_width, high_width, low_R, high_R, low_hc, high_hc, Z_TD, N = 10000) {
  
  # Sample width and recharge uniformly across measured observation windows
  width_samples <- runif(N, low_width, high_width)
  
  r_min <- min(low_R, high_R)
  r_max <- max(low_R, high_R)
  R_samples <- runif(N, r_min, r_max)
  
  # Log-Uniform Hydraulic Conductivity sampling to evaluate scale orders evenly
  hc_min <- min(low_hc, high_hc)
  hc_max <- max(low_hc, high_hc)
  HC_samples <- exp(runif(N, log(hc_min), log(hc_max)))
  
  Z_samples <- numeric(N)
  
  # Structural Vectorization for optimized memory throughput
  C_vec <- get_C_factor(width_samples)
  B_vec <- get_B_factor(width_samples)
  Y_vec <- get_Y_factor(width_samples)
  
  for (i in seq_len(N)) {
    k_val <- interpolate_K(width_samples[i], HC_samples[i])
    Z_samples[i] <- Bailey_Model(
      Y = Y_vec[i], Z_TD = Z_TD, R = R_samples[i], 
      B = B_vec[i], K = k_val, C = C_vec[i]
    )
  }
  
  return(data.frame(width = width_samples, R = R_samples, Z_TD = rep(Z_TD, N),
                    HC = HC_samples, Z_MAX = Z_samples))
}

# Sub-Iteration Vector Parser
MC_Bailey_island <- function(island_row, N = 10000) {
  df <- MC_Bailey(
    low_width  = as.numeric(island_row$Width.1),
    high_width = as.numeric(island_row$Width.4),
    low_R      = as.numeric(island_row$CDR_R),
    high_R     = as.numeric(island_row$NASA_R),
    low_hc     = 25,  # Upper Boundary
    high_hc    = 400, # Lower Boundary
    Z_TD       = as.numeric(island_row$Z_TD),
    N          = N
  )
  df$Island  <- island_row$Island
  df$Region  <- island_row$Region
  df$Z_MAX_O <- as.numeric(island_row$Z_MAX_O)
  return(df)
}


# ==============================================================================
# 4. EXECUTION PIPELINE & STATISTICAL HIGHLIGHT CALCULATIONS
# ==============================================================================

# Executing pipeline strictly across the validated variable 'island_features2'
results_all <- island_features2 %>%  
  split(.$Island) %>%   
  map(MC_Bailey_island)

results_all_df <- bind_rows(results_all)

# Save complete raw output dataset
write.csv(results_all_df, "Data_Files/Monte_Carlo_results/Monte_Carlo_Results.csv", row.names = FALSE)
results_all_df <- read.csv("Data_Files/Monte_Carlo_results/Monte_Carlo_Results.csv")

# Generate Dynamic Validation Filter Vectors to avoid typing names out twice
validation_islands <- unique(island_features2$Island)

capture_summary <- results_all_df %>%
  filter(Island %in% validation_islands) %>%
  group_by(Region, Island) %>%
  summarise(
    Observed_Z = first(Z_MAX_O),
    Median_Sim = median(Z_MAX, na.rm = TRUE),
    P2_5       = quantile(Z_MAX, 0.025, na.rm = TRUE), 
    P97_5      = quantile(Z_MAX, 0.975, na.rm = TRUE), 
    Captured_95CI = (Observed_Z >= P2_5 & Observed_Z <= P97_5),
    .groups = "drop"
  ) %>%
  # Automatically arrange islands logically by thickness within facet panels
  mutate(Island = reorder(Island, Median_Sim))


# Balanced left/right column assignment so both sides are similar height
region_column <- c(
  "Maldives"       = "A", "Cook_Islands"  = "A", "North_Pacific" = "A",
  "Kiribati"       = "B", "Marshalls"     = "B", "Micronesia"    = "B",
  "East_Indian"    = "B", "Central_Indian"= "B"
)

capture_summary <- capture_summary %>%
  mutate(Column = region_column[as.character(Region)])

# Fixed color mapping so regions stay visually consistent across both panels
region_levels    <- c("Central_Indian", "Cook_Islands", "East_Indian", "Kiribati",
                      "Maldives", "Marshalls", "Micronesia", "North_Pacific")
region_mc_colors <- setNames(scales::hue_pal()(length(region_levels)), region_levels)

# Log calculation metrics to check the pipeline status
print(capture_summary)
write.csv(capture_summary, "Data_Files/Monte_Carlo_results/capture_summary.csv", row.names = FALSE)
capture_summary <- read.csv("Data_Files/Monte_Carlo_results/capture_summary.csv")
# ==============================================================================
# 5. JOURNAL-QUALITY GRAPHICS GENERATION — TWO-COLUMN LAYOUT
# ==============================================================================

# Shorten region labels so strip boxes shrink automatically
capture_summary <- capture_summary %>%
  mutate(Region = recode(Region,
                         "Central_Indian" = "C. Indian",
                         "Cook_Islands"   = "Cook Is.",
                         "East_Indian"    = "E. Indian",
                         "North_Pacific"  = "N. Pacific"
                         # Kiribati, Maldives, Marshalls, Micronesia stay as-is — already short
  ))

# Update palette names to match the shortened labels
names(region_mc_colors) <- recode(names(region_mc_colors),
                                  "Central_Indian" = "C. Indian",
                                  "Cook_Islands"   = "Cook Is.",
                                  "East_Indian"    = "E. Indian",
                                  "North_Pacific"  = "N. Pacific"
)

make_validation_panel <- function(data) {
  ggplot(data, aes(y = Island)) +
    geom_errorbar(aes(xmin = P2_5, xmax = P97_5, color = Region),
                  width = 0.50, linewidth = 1, alpha = 0.95, orientation = "y") +
    geom_point(aes(x = Median_Sim), color = "black", size = 1.6, shape = 1, stroke = 0.9) +
    geom_point(aes(x = Observed_Z), color = "red", size = 2.3, shape = 18) +
    facet_grid(Region ~ ., scales = "free_y", space = "free_y") +
    scale_color_manual(values = region_mc_colors) +
    xlim(0, 30) +
    theme_linedraw(base_size = 9) +
    labs(x = "Freshwater Lens Thickness Z\u2098\u2090\u2093 (m)", y = NULL) +
    theme(
      legend.position    = "none",
      strip.text.y       = element_text(angle = 0, face = "bold", size = 6, color = "grey20",
                                        margin = margin(l = 2, r = 2)),
      strip.background   = element_rect(fill = "grey90", color = "grey70"),
      axis.text.y        = element_text(size = 6.5, color = "black", margin = margin(r = 2)),
      axis.ticks.length.y = unit(1, "pt"),
      axis.text.x        = element_text(size = 8, color = "black"),
      axis.title.x       = element_text(face = "bold", size = 9, margin = margin(t = 6)),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(color = "grey92", linetype = "dashed"),
      panel.spacing.y    = unit(2, "pt"),
      plot.margin        = margin(2, 2, 2, 2)
    )
}

panel_A <- make_validation_panel(capture_summary %>% filter(Column == "A"))
panel_B <- make_validation_panel(capture_summary %>% filter(Column == "B"))

validation_plot <- panel_A + panel_B +
  plot_layout(ncol = 2) +
  plot_annotation(
    title    = "Freshwater Lens Model Calibration & Validation Matrix",
    subtitle = "Bars = 95% simulated interval | Circles = model median | Red diamonds = field observation",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 8, color = "grey30")
    )
  )

ggsave(
  filename = "R_Plots/Paper_MC_Validation.png",
  plot     = validation_plot,
  width    = 15,     # ERL double-column width
  height   = 18,     # tall enough for ~25 rows/column at readable size
  units    = "cm",
  dpi      = 600
)

ggsave(
  filename = "R_Plots/Paper_MC_Validation.pdf",
  plot     = validation_plot,
  width    = 15,
  height   = 18,
  units    = "cm"
)


print(validation_plot)
