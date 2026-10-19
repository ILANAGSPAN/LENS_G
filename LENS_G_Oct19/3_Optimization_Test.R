

# ==============================================================================
# BAILEY ET AL. 2008 - FRESHWATER LENS THICKNESS MODEL
# Full optimization script:
#   Model A = Strict Bailey
#   Model B = Global recharge optimization
#   Model C = Regional recharge optimization
#   Model D = Regional recharge + regional K scalar
#
# PRECIPITATION INPUT:
#   MSWEP_R, JAXA_R, NASA_R, CDR_R, and ERA5_R in FOR_PROJECTION.csv
#   are product-specific recharge estimates in m/yr, calculated as:
#
#       Recharge = annual precipitation × 0.5
#
#   These columns are therefore converted back to precipitation by
#   multiplying by 2, and the five reconstructed precipitation estimates
#   are averaged to create the model precipitation input:
#
#       Ensemble.Precip =
#           mean(MSWEP_R, JAXA_R, NASA_R, CDR_R, ERA5_R) × 2
#
#   Model A then applies the Bailey recharge fraction of 0.5.
#   Models B-D optimize recharge fractions applied to Ensemble.Precip.
#
# MODEL D:
#   K_effective = K_orientation × regional K scalar
#
#   K_orientation:
#       windward = 400 m/day
#       leeward  = 50 m/day
#       offset   = 200 m/day
#
# ==============================================================================


# ==============================================================================
# PACKAGES
# ==============================================================================

library(dplyr)
library(ggplot2)
library(purrr)
library(minpack.lm)
library(DEoptim)
library(DescTools)
library(patchwork)
library(cowplot)
library(tidyr)

setwd("Q:/Hydro/projects/indopacom_PRIVATE/ILANA")
setwd("C:/Users/Owner/OneDrive/Desktop/Github") 

load("C:/Users/Owner/OneDrive/Desktop/MS_Package_Oct7/Oct7.RData")

# ==============================================================================
# SECTION 1: DIGITIZED FACTOR CURVES
# ==============================================================================

df_b <- data.frame(
  island_width = c(0, 25, 70, 120, 190, 250, 290, 400, 500, 700,
                   800, 900, 950, 1000, 1050, 1100, 1150, 1200,
                   1500, 1800, 2000, 2200, 2300, 2500, 2600),
  b_factor     = c(6.9, 5.2, 4, 3, 2, 1.3, 1, 0.6, 0.3, 0.25,
                   0.15, 0.1, 0.08, 0.06, 0.045, 0.032, 0.021, 0.015,
                   0.010, 0.006, 0.004, 0.003, 0.001, 0.0001, 0.0001)
)

df_c <- data.frame(
  island_width = c(0, 125, 250, 275, 500, 750, 800, 900, 1000, 1050,
                   1200, 1500, 1600, 2000, 2200, 2500, 2600),
  c_factor     = c(0.853, 0.885, 0.91, 0.92, 0.955, 0.98, 0.985, 0.99,
                   0.995, 0.996, 1, 1, 1, 1, 1, 1, 1)
)

df_y <- data.frame(
  y_factor_widths = c(0, 500, 1000, 1500, 2000, 2500),
  y_factor        = c(-1, -9.9, -14.9, -17.8, -19.9, -21.2)
)

island_width_k <- data.frame(
  width = c(150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750,
            800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1300, 1400,
            1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500),
  
  k25  = c(1.6, 1.6, 1.58, 1.54, 1.47, 1.36, 1.24, 1.13, 1.07, 1.03, 1.02,
           1.01, 1.01, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1),
  
  k50  = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
  
  k75  = c(0.74, 0.75, 0.77, 0.79, 0.81, 0.83, 0.86, 0.88, 0.91, 0.93, 0.94,
           0.96, 0.97, 0.98, 0.98, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99,
           0.99, 0.99, 0.99, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
  
  k100 = c(0.61, 0.62, 0.63, 0.66, 0.68, 0.72, 0.75, 0.79, 0.83, 0.87, 0.9,
           0.92, 0.94, 0.95, 0.96, 0.97, 0.98, 0.98, 0.98, 0.98, 0.99, 0.99,
           0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99,
           0.99, 0.99),
  
  k125 = c(0.52, 0.53, 0.55, 0.57, 0.6, 0.63, 0.67, 0.71, 0.76, 0.8, 0.84,
           0.87, 0.9, 0.92, 0.94, 0.95, 0.96, 0.97, 0.97, 0.97, 0.98, 0.98,
           0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98,
           0.98, 0.98),
  
  k150 = c(0.46, 0.47, 0.49, 0.51, 0.53, 0.57, 0.6, 0.65, 0.69, 0.74, 0.78,
           0.82, 0.86, 0.88, 0.91, 0.93, 0.94, 0.95, 0.96, 0.96, 0.97, 0.97,
           0.97, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98,
           0.98, 0.98),
  
  k175 = c(0.41, 0.42, 0.44, 0.46, 0.48, 0.51, 0.55, 0.59, 0.63, 0.68, 0.72,
           0.77, 0.81, 0.84, 0.87, 0.89, 0.91, 0.93, 0.94, 0.95, 0.96, 0.96,
           0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97,
           0.97, 0.97),
  
  k200 = c(0.38, 0.39, 0.4, 0.42, 0.44, 0.47, 0.5, 0.53, 0.57, 0.62, 0.66,
           0.71, 0.75, 0.79, 0.83, 0.86, 0.88, 0.9, 0.92, 0.93, 0.94, 0.95,
           0.96, 0.96, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97,
           0.97, 0.97),
  
  k225 = c(0.35, 0.36, 0.37, 0.39, 0.4, 0.43, 0.45, 0.49, 0.53, 0.57, 0.61,
           0.66, 0.7, 0.74, 0.78, 0.82, 0.85, 0.87, 0.89, 0.91, 0.92, 0.93,
           0.94, 0.95, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96,
           0.96, 0.96),
  
  k250 = c(0.33, 0.33, 0.34, 0.36, 0.37, 0.4, 0.42, 0.45, 0.48, 0.52, 0.56,
           0.61, 0.65, 0.7, 0.74, 0.77, 0.81, 0.84, 0.86, 0.88, 0.9, 0.91,
           0.93, 0.94, 0.95, 0.95, 0.95, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96,
           0.96, 0.96),
  
  k275 = c(0.31, 0.31, 0.32, 0.34, 0.35, 0.37, 0.39, 0.42, 0.45, 0.48, 0.52,
           0.56, 0.61, 0.65, 0.69, 0.73, 0.77, 0.8, 0.83, 0.85, 0.87, 0.89,
           0.91, 0.93, 0.94, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.96,
           0.96, 0.96),
  
  k300 = c(0.29, 0.3, 0.31, 0.32, 0.33, 0.35, 0.37, 0.39, 0.42, 0.45, 0.49,
           0.52, 0.56, 0.61, 0.65, 0.69, 0.73, 0.76, 0.79, 0.82, 0.85, 0.87,
           0.9, 0.92, 0.94, 0.94, 0.95, 0.95, 0.96, 0.96, 0.96, 0.96, 0.96,
           0.96, 0.96),
  
  k325 = c(0.27, 0.28, 0.29, 0.3, 0.31, 0.33, 0.35, 0.37, 0.4, 0.42, 0.46,
           0.49, 0.53, 0.57, 0.61, 0.65, 0.69, 0.73, 0.76, 0.79, 0.82, 0.85,
           0.89, 0.92, 0.94, 0.95, 0.96, 0.97, 0.97, 0.97, 0.98, 0.98, 0.98,
           0.98, 0.98),
  
  k400 = c(0.24, 0.25, 0.25, 0.26, 0.27, 0.29, 0.3, 0.32, 0.34, 0.37, 0.39,
           0.42, 0.45, 0.49, 0.53, 0.56, 0.6, 0.64, 0.68, 0.71, 0.74, 0.77,
           0.82, 0.86, 0.89, 0.91, 0.92, 0.93, 0.94, 0.94, 0.95, 0.95, 0.95,
           0.95, 0.95)
)

K_VALUES <- c(
  25, 50, 75, 100, 125, 150, 175, 200,
  225, 250, 275, 300, 325, 400
)


# ==============================================================================
# SECTION 2: SHAPE FACTOR FUNCTIONS
# ==============================================================================

get_C <- function(width) {
  1 - 0.15 * exp(-0.00288 * width)
}

get_B <- function(width) {
  5 * exp(-0.00461 * width)
}

get_Y <- function(width) {
  -20 * (1 - exp(-0.001386 * width))
}


# ==============================================================================
# SECTION 3: K CORRECTION FUNCTION
# ==============================================================================

get_K_correction <- function(width, K_mday) {
  
  if (is.na(width) || is.na(K_mday) || K_mday <= 0) {
    return(NA_real_)
  }
  
  width_clamped <- pmin(pmax(width, 150), 2500)
  K_clamped     <- pmin(pmax(K_mday, 25), 400)
  
  corrections <- sapply(
    paste0("k", K_VALUES),
    function(col) {
      approx(
        island_width_k$width,
        island_width_k[[col]],
        xout = width_clamped,
        rule = 2
      )$y
    }
  )
  
  approx(
    K_VALUES,
    corrections,
    xout = K_clamped,
    rule = 2
  )$y
}


# ==============================================================================
# SECTION 4: ORIENTATION → BASELINE K
# ==============================================================================

get_K_orientation <- function(orient) {
  
  case_when(
    !is.na(orient) & orient == "windward" ~ 400,
    !is.na(orient) & orient == "leeward"  ~ 50,
    !is.na(orient) & orient == "offset"   ~ 200,
    TRUE                                  ~ 200
  )
}


# ==============================================================================
# SECTION 5: CORE MODEL FUNCTION
#
# rainfall_m    = annual precipitation in m/yr
# recharge_frac = fraction of precipitation becoming recharge
# K_mday        = physical hydraulic conductivity in m/day
#
# The recharge fraction is applied HERE.
#
# Therefore:
#   Model A: recharge_frac = 0.5
#   Model B: one global optimized recharge fraction
#   Model C: one optimized fraction per region
#   Model D: one optimized fraction per region
#
# ==============================================================================

bailey_core <- function(
    rainfall_m,
    width,
    Z_TD,
    recharge_frac,
    K_mday
) {
  
  if (any(is.na(c(
    rainfall_m,
    width,
    Z_TD,
    recharge_frac
  )))) {
    return(NA_real_)
  }
  
  if (is.na(K_mday) || K_mday <= 0) {
    return(NA_real_)
  }
  
  # Annual recharge (m/yr)
  R <- rainfall_m * recharge_frac
  
  B <- get_B(width)
  C <- get_C(width)
  Y <- get_Y(width)
  
  K_corr <- get_K_correction(width, K_mday)
  
  if (is.na(K_corr)) {
    return(NA_real_)
  }
  
  Z_MAX <- (
    ((Z_TD - Y) * R) / (B + R) + Y
  ) * K_corr * C
  
  as.numeric(Z_MAX)
}


# ==============================================================================
# SECTION 6: DATA LOADING AND ENSEMBLE PRECIPITATION
# ==============================================================================

island_features <- read.csv(
  "Data_Files/FOR_PROJECTION.csv"
)

# Convert blank strings to NA
island_features <- island_features %>%
  mutate(
    across(
      everything(),
      ~ na_if(as.character(.), "")
    )
  )

# Convert numeric columns back to numeric where appropriate
numeric_columns <- c(
  "MSWEP_R",
  "JAXA_R",
  "NASA_R",
  "CDR_R",
  "ERA5_R",
  "Avg.Width",
  "Z_TD",
  "Z_MAX_O",
  "Z_MAX_B",
  "Z_MAX_P",
  "Avg.Precip"
)

numeric_columns <- intersect(
  numeric_columns,
  names(island_features)
)

island_features <- island_features %>%
  mutate(
    across(
      all_of(numeric_columns),
      as.numeric
    )
  )


# ------------------------------------------------------------------------------
# IMPORTANT:
#
# The *_R columns are already recharge estimates:
#
#   *_R = annual precipitation × 0.5
#
# Therefore:
#
#   precipitation = *_R × 2
#
# We reconstruct each product's precipitation and then average the five
# precipitation estimates.
# ------------------------------------------------------------------------------

required_R_columns <- c(
  "MSWEP_R",
  "JAXA_R",
  "NASA_R",
  "CDR_R",
  "ERA5_R"
)

missing_R_columns <- setdiff(
  required_R_columns,
  names(island_features)
)

if (length(missing_R_columns) > 0) {
  stop(
    paste(
      "The following recharge columns are missing from FOR_PROJECTION.csv:",
      paste(missing_R_columns, collapse = ", ")
    )
  )
}


# Reconstruct precipitation for each product and calculate the ensemble.
# If one product is unavailable, calculate the ensemble from the
# remaining products, provided that at least four of the five products
# are available.

island_features <- island_features %>%
  mutate(
    
    MSWEP_Precip = MSWEP_R * 2,
    JAXA_Precip  = JAXA_R  * 2,
    NASA_Precip  = NASA_R  * 2,
    CDR_Precip   = CDR_R   * 2,
    ERA5_Precip  = ERA5_R  * 2,
    
    n_precip_products = rowSums(
      !is.na(
        cbind(
          MSWEP_Precip,
          JAXA_Precip,
          NASA_Precip,
          CDR_Precip,
          ERA5_Precip
        )
      )
    ),
    
    Ensemble.Precip = ifelse(
      n_precip_products >= 4,
      rowMeans(
        cbind(
          MSWEP_Precip,
          JAXA_Precip,
          NASA_Precip,
          CDR_Precip,
          ERA5_Precip
        ),
        na.rm = TRUE
      ),
      NA_real_
    )
  )

# ------------------------------------------------------------------------------
# Diagnostic check
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("ENSEMBLE PRECIPITATION CHECK\n")
cat("============================================================\n")

print(
  island_features %>%
    dplyr::select(
      Island.ID,
      Region,
      Island,
      MSWEP_R,
      JAXA_R,
      NASA_R,
      CDR_R,
      ERA5_R,
      Ensemble.Precip
    ) %>%
    head()
)


# Check that ensemble precipitation is plausible
ensemble_summary <- island_features %>%
  summarise(
    n = sum(is.finite(Ensemble.Precip)),
    mean = mean(Ensemble.Precip, na.rm = TRUE),
    min = min(Ensemble.Precip, na.rm = TRUE),
    max = max(Ensemble.Precip, na.rm = TRUE),
    missing = sum(!is.finite(Ensemble.Precip))
  )

print(ensemble_summary)


# ------------------------------------------------------------------------------
# Filter regions
# ------------------------------------------------------------------------------

final_results <- island_features %>%
  mutate(
    Region = gsub(
      " ",
      "_",
      trimws(Region)
    )
  ) %>%
  filter(
    Region %in% c(
      "Marshalls",
      "FSM",
      "Maldives",
      "Central_Indian",
      "East_Indian",
      "Kiribati",
      "North_Pacific",
      "Cook_Islands"
    )
  )


# ------------------------------------------------------------------------------
# Diagnostic: missing key inputs
# ------------------------------------------------------------------------------

missing_check <- final_results %>%
  mutate(row = seq_len(nrow(final_results))) %>%
  filter(
    is.na(Avg.Width) |
      is.na(Ensemble.Precip) |
      is.na(Z_TD) |
      is.na(Orientation)
  ) %>%
  dplyr::select(
    row,
    Region,
    Island,
    Avg.Width,
    Ensemble.Precip,
    Orientation,
    Z_TD
  )

if (nrow(missing_check) > 0) {
  
  cat(
    "\nWARNING: The following islands have missing inputs ",
    "and will return NA predictions:\n"
  )
  
  print(missing_check)
}


regions_present <- unique(final_results$Region)

n_regions <- length(regions_present)


# ==============================================================================
# SECTION 7: REGION PLOTTING INFORMATION
# ==============================================================================

region_colors <- c(
  "Marshalls"      = "firebrick",
  "FSM"            = "cyan2",
  "Maldives"       = "darkorchid",
  "Central_Indian" = "green",
  "East_Indian"    = "orange",
  "Kiribati"       = "darkcyan",
  "North_Pacific"  = "blue",
  "Cook_Islands"   = "magenta"
)

region_shapes <- c(
  "Marshalls"      = 16,
  "FSM"            = 17,
  "Maldives"       = 15,
  "Central_Indian" = 18,
  "East_Indian"    = 1,
  "Kiribati"       = 2,
  "North_Pacific"  = 0,
  "Cook_Islands"   = 5
)

region_labels <- c(
  "Marshalls"      = "Marshall Islands",
  "FSM"            = "Micronesia",
  "Maldives"       = "Maldives",
  "Central_Indian" = "Central Indian Ocean",
  "East_Indian"    = "East Indian Ocean",
  "Kiribati"       = "Kiribati",
  "North_Pacific"  = "North Pacific",
  "Cook_Islands"   = "Cook Islands"
)


# ==============================================================================
# SECTION 8: MODEL A — STRICT BAILEY BASELINE
#
# Bailey assumption:
#   recharge = 50% of annual precipitation
#
# K is determined solely by orientation.
# ==============================================================================

final_results <- final_results %>%
  mutate(
    
    Z_MAX_modelA = pmap_dbl(
      list(
        Ensemble.Precip,
        Avg.Width,
        Z_TD,
        Orientation
      ),
      
      function(
    rain,
    width,
    ztd,
    orient
      ) {
        
        K_mday <- get_K_orientation(orient)
        
        bailey_core(
          rain,
          width,
          ztd,
          recharge_frac = 0.5,
          K_mday = K_mday
        )
      }
    )
  )


metrics_A <- final_results %>%
  filter(
    is.finite(Z_MAX_modelA) &
      is.finite(Z_MAX_O)
  ) %>%
  summarise(
    Model = "A: Strict Bailey",
    
    RMSE = sqrt(
      mean(
        (Z_MAX_modelA - Z_MAX_O)^2
      )
    ),
    
    MAE = mean(
      abs(
        Z_MAX_modelA - Z_MAX_O
      )
    ),
    
    Bias = mean(
      Z_MAX_modelA - Z_MAX_O
    ),
    
    R2 = cor(
      Z_MAX_O,
      Z_MAX_modelA
    )^2,
    
    CCC = DescTools::CCC(
      Z_MAX_O,
      Z_MAX_modelA
    )rho.cest
  )

print(metrics_A)


# ==============================================================================
# SECTION 9: MODEL B — GLOBAL RECHARGE OPTIMIZATION
# ==============================================================================

obj_fn_B <- function(
    params,
    data
) {
  
  rf <- params[1]
  
  pred <- pmap_dbl(
    list(
      data$Ensemble.Precip,
      data$Avg.Width,
      data$Z_TD,
      data$Orientation
    ),
    
    function(
    rain,
    width,
    ztd,
    orient
    ) {
      
      K_mday <- get_K_orientation(orient)
      
      bailey_core(
        rain,
        width,
        ztd,
        rf,
        K_mday
      )
    }
  )
  
  valid <- is.finite(pred) &
    is.finite(data$Z_MAX_O)
  
  if (sum(valid) < 5) {
    return(1e10)
  }
  
  sqrt(
    mean(
      (pred[valid] - data$Z_MAX_O[valid])^2
    )
  )
}


cat("\nRunning Model B optimization...\n")

library(DEoptim)

fit_B <- DEoptim(
  fn = obj_fn_B,
  
  lower = 0.1,
  upper = 0.9,
  
  data = final_results,
  
  control = DEoptim.control(
    itermax = 300,
    NP = 30,
    trace = 100
  )
)

rf_global <- fit_Boptimbestmem

cat(
  "\nModel B — Optimal global recharge fraction:",
  round(rf_global, 3),
  "\n"
)


final_results <- final_results %>%
  mutate(
    
    Z_MAX_modelB = pmap_dbl(
      list(
        Ensemble.Precip,
        Avg.Width,
        Z_TD,
        Orientation
      ),
      
      function(
    rain,
    width,
    ztd,
    orient
      ) {
        
        K_mday <- get_K_orientation(orient)
        
        bailey_core(
          rain,
          width,
          ztd,
          rf_global,
          K_mday
        )
      }
    )
  )


metrics_B <- final_results %>%
  filter(
    is.finite(Z_MAX_modelB) &
      is.finite(Z_MAX_O)
  ) %>%
  summarise(
    
    Model = "B: Global Recharge",
    
    RMSE = sqrt(
      mean(
        (Z_MAX_modelB - Z_MAX_O)^2
      )
    ),
    
    MAE = mean(
      abs(
        Z_MAX_modelB - Z_MAX_O
      )
    ),
    
    Bias = mean(
      Z_MAX_modelB - Z_MAX_O
    ),
    
    R2 = cor(
      Z_MAX_O,
      Z_MAX_modelB
    )^2,
    
    CCC = DescTools::CCC(
      Z_MAX_O,
      Z_MAX_modelB
    )rho.cest
  )

print(metrics_B)


# ==============================================================================
# SECTION 10: MODEL C — REGIONAL RECHARGE OPTIMIZATION
# ==============================================================================

obj_fn_C <- function(
    params,
    data
) {
  
  rf_by_region <- setNames(
    params,
    paste0(
      "rf_",
      regions_present
    )
  )
  
  pred <- pmap_dbl(
    list(
      data$Ensemble.Precip,
      data$Avg.Width,
      data$Z_TD,
      data$Orientation,
      data$Region
    ),
    
    function(
    rain,
    width,
    ztd,
    orient,
    region
    ) {
      
      rf <- rf_by_region[
        paste0(
          "rf_",
          region
        )
      ]
      
      K_mday <- get_K_orientation(
        orient
      )
      
      bailey_core(
        rain,
        width,
        ztd,
        rf,
        K_mday
      )
    }
  )
  
  valid <- is.finite(pred) &
    is.finite(data$Z_MAX_O)
  
  if (sum(valid) < 10) {
    return(1e10)
  }
  
  sqrt(
    mean(
      (pred[valid] - data$Z_MAX_O[valid])^2
    )
  )
}


cat("\nRunning Model C optimization...\n")

fit_C <- DEoptim(
  fn = obj_fn_C,
  
  lower = rep(
    0.1,
    n_regions
  ),
  
  upper = rep(
    0.9,
    n_regions
  ),
  
  data = final_results,
  
  control = DEoptim.control(
    itermax = 500,
    NP = 10 * n_regions,
    trace = 100
  )
)


rf_regional <- setNames(
  fit_Coptimbestmem,
  paste0(
    "rf_",
    regions_present
  )
)


cat(
  "\nModel C — Optimal recharge fractions by region:\n"
)

print(
  round(
    rf_regional,
    3
  )
)


final_results <- final_results %>%
  mutate(
    
    Z_MAX_modelC = pmap_dbl(
      list(
        Ensemble.Precip,
        Avg.Width,
        Z_TD,
        Orientation,
        Region
      ),
      
      function(
    rain,
    width,
    ztd,
    orient,
    region
      ) {
        
        rf <- rf_regional[
          paste0(
            "rf_",
            region
          )
        ]
        
        K_mday <- get_K_orientation(
          orient
        )
        
        bailey_core(
          rain,
          width,
          ztd,
          rf,
          K_mday
        )
      }
    )
  )


metrics_C <- final_results %>%
  filter(
    is.finite(Z_MAX_modelC) &
      is.finite(Z_MAX_O)
  ) %>%
  summarise(
    
    Model = "C: Regional Recharge",
    
    RMSE = sqrt(
      mean(
        (Z_MAX_modelC - Z_MAX_O)^2
      )
    ),
    
    MAE = mean(
      abs(
        Z_MAX_modelC - Z_MAX_O
      )
    ),
    
    Bias = mean(
      Z_MAX_modelC - Z_MAX_O
    ),
    
    R2 = cor(
      Z_MAX_O,
      Z_MAX_modelC
    )^2,
    
    CCC = DescTools::CCC(
      Z_MAX_O,
      Z_MAX_modelC
    )rho.cest
  )

print(metrics_C)


# ==============================================================================
# SECTION 11: MODEL D — REGIONAL RECHARGE + REGIONAL K SCALAR
#
# K_effective =
#     K_orientation × regional_K_scalar
#
# Orientation values remain:
#     windward = 400 m/day
#     leeward  = 50 m/day
#     offset   = 200 m/day
#
# ==============================================================================

obj_fn_D <- function(
    params,
    data
) {
  
  rf_by_region <- setNames(
    params[1:n_regions],
    paste0(
      "rf_",
      regions_present
    )
  )
  
  Kscale_by_region <- setNames(
    params[
      (n_regions + 1):(2 * n_regions)
    ],
    paste0(
      "Ks_",
      regions_present
    )
  )
  
  pred <- pmap_dbl(
    list(
      data$Ensemble.Precip,
      data$Avg.Width,
      data$Z_TD,
      data$Orientation,
      data$Region
    ),
    
    function(
    rain,
    width,
    ztd,
    orient,
    region
    ) {
      
      rf <- rf_by_region[
        paste0(
          "rf_",
          region
        )
      ]
      
      K_scalar <- Kscale_by_region[
        paste0(
          "Ks_",
          region
        )
      ]
      
      K_mday <- get_K_orientation(
        orient
      ) * K_scalar
      
      K_mday <- pmin(
        pmax(
          K_mday,
          25
        ),
        400
      )
      
      bailey_core(
        rain,
        width,
        ztd,
        rf,
        K_mday
      )
    }
  )
  
  valid <- is.finite(pred) &
    is.finite(data$Z_MAX_O)
  
  if (sum(valid) < 10) {
    return(1e10)
  }
  
  sqrt(
    mean(
      (pred[valid] - data$Z_MAX_O[valid])^2
    )
  )
}


cat("\nRunning Model D optimization...\n")


fit_D <- DEoptim(
  fn = obj_fn_D,
  
  lower = c(
    rep(
      0.1,
      n_regions
    ),
    
    rep(
      0.25,
      n_regions
    )
  ),
  
  upper = c(
    rep(
      0.9,
      n_regions
    ),
    
    rep(
      4.0,
      n_regions
    )
  ),
  
  data = final_results,
  
  control = DEoptim.control(
    itermax = 600,
    NP = 10 * (2 * n_regions),
    trace = 100
  )
)


params_D <- fit_Doptimbestmem


rf_D <- setNames(
  params_D[
    1:n_regions
  ],
  paste0(
    "rf_",
    regions_present
  )
)


Kscale_D <- setNames(
  params_D[
    (n_regions + 1):(2 * n_regions)
  ],
  paste0(
    "Ks_",
    regions_present
  )
)


cat(
  "\nModel D — Optimal recharge fractions:\n"
)

print(
  round(
    rf_D,
    3
  )
)


cat(
  "\nModel D — Regional K scalars:\n"
)

print(
  round(
    Kscale_D,
    3
  )
)


# ------------------------------------------------------------------------------
# Effective K table
# ------------------------------------------------------------------------------

K_effective_table <- data.frame(
  
  Region = regions_present,
  
  Scalar = round(
    as.numeric(
      Kscale_D
    ),
    3
  ),
  
  K_windward_effective = round(
    pmin(
      pmax(
        400 *
          as.numeric(Kscale_D),
        25
      ),
      400
    ),
    1
  ),
  
  K_leeward_effective = round(
    pmin(
      pmax(
        50 *
          as.numeric(Kscale_D),
        25
      ),
      400
    ),
    1
  ),
  
  K_offset_effective = round(
    pmin(
      pmax(
        200 *
          as.numeric(Kscale_D),
        25
      ),
      400
    ),
    1
  )
)

print(
  K_effective_table
)

write.csv(
  K_effective_table,
  "Data_Files/ModelD_effective_K.csv",
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# Model D predictions
# ------------------------------------------------------------------------------

final_results <- final_results %>%
  mutate(
    
    Z_MAX_modelD = pmap_dbl(
      list(
        Ensemble.Precip,
        Avg.Width,
        Z_TD,
        Orientation,
        Region
      ),
      
      function(
    rain,
    width,
    ztd,
    orient,
    region
      ) {
        
        rf <- rf_D[
          paste0(
            "rf_",
            region
          )
        ]
        
        K_scalar <- Kscale_D[
          paste0(
            "Ks_",
            region
          )
        ]
        
        K_mday <- get_K_orientation(
          orient
        ) * K_scalar
        
        K_mday <- pmin(
          pmax(
            K_mday,
            25
          ),
          400
        )
        
        bailey_core(
          rain,
          width,
          ztd,
          rf,
          K_mday
        )
      }
    )
  )


metrics_D <- final_results %>%
  filter(
    is.finite(Z_MAX_modelD) &
      is.finite(Z_MAX_O)
  ) %>%
  summarise(
    
    Model = "D: Regional Recharge + K Scalar",
    
    RMSE = sqrt(
      mean(
        (Z_MAX_modelD - Z_MAX_O)^2
      )
    ),
    
    MAE = mean(
      abs(
        Z_MAX_modelD - Z_MAX_O
      )
    ),
    
    Bias = mean(
      Z_MAX_modelD - Z_MAX_O
    ),
    
    R2 = cor(
      Z_MAX_O,
      Z_MAX_modelD
    )^2,
    
    CCC = DescTools::CCC(
      Z_MAX_O,
      Z_MAX_modelD
    )rho.cest
  )

print(metrics_D)


# ==============================================================================
# SECTION 12: MODEL COMPARISON TABLE
# ==============================================================================

all_metrics <- bind_rows(
  metrics_A,
  metrics_B,
  metrics_C,
  metrics_D
) %>%
  mutate(
    across(
      where(is.numeric),
      ~ round(., 3)
    )
  )

print(all_metrics)

write.csv(
  all_metrics,
  "Data_Files/model_comparison_metrics.csv",
  row.names = FALSE
)


# ==============================================================================
# SECTION 13: REGIONAL PARAMETER TABLE
# ==============================================================================

parameter_plot_df <- tibble(
  
  Region = regions_present,
  
  Recharge = as.numeric(
    rf_D
  ),
  
  K_scalar = as.numeric(
    Kscale_D
  )
  
) %>%
  mutate(
    Region_label = region_labels[Region]
  ) %>%
  arrange(
    Recharge
  ) %>%
  mutate(
    Region_label = factor(
      Region_label,
      levels = Region_label
    )
  )

# ==============================================================================
# SECTION 8: VISUALIZATION — MODEL PROGRESSION
# ==============================================================================

library(ggplot2)

model_progression <- data.frame(
  
  Model = c(
    "A: Strict Bailey\n(RF=0.5, K=orientation)",
    "B: Global RF\n(RF opt, K=orientation)",
    "C: Regional RF\n(RF opt, K=orientation)",
    "D: Regional RF + K\n(Both opt)"
  ),
  
  RMSE = c(
    metrics_A$RMSE,
    metrics_B$RMSE,
    metrics_C$RMSE,
    metrics_D$RMSE
  ),
  
  R2 = c(
    metrics_A$R2,
    metrics_B$R2,
    metrics_C$R2,
    metrics_D$R2
  ),
  
  CCC = c(
    metrics_A$CCC,
    metrics_B$CCC,
    metrics_C$CCC,
    metrics_D$CCC
  )
)

# Reorder for plot
model_progression$Model <- factor(
  model_progression$Model,
  levels = model_progression$Model
)

p_rmse <- ggplot(
  model_progression,
  aes(
    x = Model,
    y = RMSE,
    fill = Model
  )
) +
  geom_col(
    color = "black",
    linewidth = 1
  ) +
  labs(
    y = "RMSE (m)",
    title = "Model Progression: RMSE"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 9
    )
  )

p_r2 <- ggplot(
  model_progression,
  aes(
    x = Model,
    y = R2,
    fill = Model
  )
) +
  geom_col(
    color = "black",
    linewidth = 1
  ) +
  labs(
    y = "R²",
    title = "Model Progression: R²"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 9
    )
  ) +
  ylim(0, 1)

p_ccc <- ggplot(
  model_progression,
  aes(
    x = Model,
    y = CCC,
    fill = Model
  )
) +
  geom_col(
    color = "black",
    linewidth = 1
  ) +
  labs(
    y = "CCC",
    title = "Model Progression: Concordance"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 9
    )
  ) +
  ylim(0, 1)

p_combined <- p_rmse / p_r2 / p_ccc

ggsave(
  "Plots/Model_Progression_ABCD.png",
  p_combined,
  width = 10,
  height = 12,
  dpi = 300
)

cat("\n✓ Saved: Model_Progression_ABCD.png\n")


# =============================================================================
# FOUR-PANEL OBSERVED-VERSUS-PREDICTED PLOT (MODELS A-D)
# Uses predictions and metrics already computed above; does not rerun model fitting.
# =============================================================================

region_order <- c("Marshalls", "FSM", "Maldives", "Central_Indian",
                  "East_Indian", "Kiribati", "North_Pacific", "Cook_Islands")
final_resultsRegion<-factor(finalresultsRegion, levels = region_order)

plot_model <- function(data, pred_col, title_str, metrics_row,
                       show_x_title = TRUE, show_y_title = TRUE) {
  
  subtitle <- paste0("RMSE = ", round(metrics_row$RMSE, 2),
                     "   R\u00b2 = ", round(metrics_row$R2, 2))
  
  p <- ggplot(data, aes(x = Z_MAX_O, y = .data[[pred_col]],
                        color = Region, shape = Region)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "red", linewidth = 0.5) +
    geom_point(size = 2.1, alpha = 0.8, stroke = 1) +
    scale_color_manual(values = region_colors,
                       labels = region_labels[levels(data$Region)],
                       name   = "Region") +
    scale_shape_manual(values = region_shapes,
                       labels = region_labels[levels(data$Region)],
                       name   = "Region") +
    coord_equal(xlim = c(0, 32), ylim = c(0, 32)) +
    labs(title = title_str, subtitle = subtitle,
         x = "Observed (m)", y = "Predicted (m)") +
    theme_bw(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 11, margin = margin(b = 2)),
      plot.subtitle = element_text(size = 9, color = "gray30"),
      plot.margin   = margin(4, 6, 4, 4),
      axis.title    = element_text(size = 10),
      axis.text     = element_text(size = 9)
    )
  
  if (!show_x_title) p <- p + theme(axis.title.x = element_blank())
  if (!show_y_title) p <- p + theme(axis.title.y = element_blank())
  
  p
}

pA <- plot_model(final_results, "Z_MAX_modelA",
                 "A: Strict Bailey", metrics_A,
                 show_x_title = FALSE, show_y_title = TRUE)
pB <- plot_model(final_results, "Z_MAX_modelB",
                 "B: Global Recharge", metrics_B,
                 show_x_title = FALSE, show_y_title = FALSE)
pC <- plot_model(final_results, "Z_MAX_modelC",
                 "C: Regional Recharge", metrics_C,
                 show_x_title = TRUE, show_y_title = TRUE)
pD <- plot_model(final_results, "Z_MAX_modelD",
                 "D: Regional + K Scalar", metrics_D,
                 show_x_title = TRUE, show_y_title = FALSE)

# patchwork combines panels AND auto-collects the one shared legend
final_four_panel <-
  (pA | pB) / (pC | pD) +
  plot_layout(guides = "collect") &
  theme(legend.position  = "bottom",
        legend.title     = element_text(size = 9, face = "bold"),
        legend.text      = element_text(size = 8),
        legend.key.size  = unit(0.35, "cm")) &
  guides(color = guide_legend(nrow = 2, byrow = TRUE,
                              override.aes = list(size = 2.5)),
         shape = guide_legend(nrow = 2, byrow = TRUE))

ggsave("C:/Users/Owner/OneDrive/Desktop/MS_Package_Oct7/R_Plots/FourPanel_AllModels.png",
       final_four_panel,
       width  = 15,
       height = 17,   # nudged down from 16 — the shorter subtitle frees vertical space
       units  = "cm",
       dpi    = 600)


ggsave("C:/Users/Owner/OneDrive/Desktop/MS_Package_Oct7/R_Plots/FourPanel_AllModels.pdf",
       final_four_panel,
       width  = 15,
       height = 17,   # nudged down from 16 — the shorter subtitle frees vertical space
       units  = "cm",
       dpi    = 600)

print(final_four_panel)



# ==============================================================================
# REGIONAL PARAMETER LOLLIPOP PLOT — MODEL D
# Uses the fitted Model D estimates above; does not rerun optimization.
# ==============================================================================

parameter_plot_df <- tibble(
  Region = regions_present,
  Recharge = as.numeric(rf_D),
  K_scalar = as.numeric(Kscale_D)) %>%
  mutate(Region_label = unname(region_labels[Region])) %>%
  arrange(Recharge) %>%
  mutate(Region_label = factor(Region_label, levels = Region_label))

p_rf <- ggplot(parameter_plot_df,
               aes(x = Recharge, y = Region_label)) +
  geom_segment(aes(x = 0, xend = Recharge,
                   y = Region_label, yend = Region_label),
               linewidth = 1) +
  geom_point(size = 3) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Recharge Fraction", x = "Recharge Fraction", y = NULL) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 9.5))

p_k <- ggplot(parameter_plot_df,
              aes(x = K_scalar, y = Region_label)) +
  geom_segment(aes(x = 0, xend = K_scalar,
                   y = Region_label, yend = Region_label),
               linewidth = 1) +
  geom_point(size = 3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_continuous(limits = c(0, 1.1), breaks = seq(0, 0.9, 0.3),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "K Scalar", x = "K Scalar (1 = Bailey baseline)", y = NULL) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 9.5))

fig_parameters <- p_rf + p_k + plot_layout(widths = c(1.3, 1))

ggsave("C:/Users/Owner/OneDrive/Desktop/MS_Package_Oct7/R_Plots/Regional_parameters.png", fig_parameters,
       width = 15, height = 8, units = "cm", dpi = 600)
ggsave("C:/Users/Owner/OneDrive/Desktop/MS_Package_Oct7/R_Plots/Regional_parameters.pdf", fig_parameters,
       width = 15, height = 8, units = "cm")
print(fig_parameters)
cat("Saved Model D regional parameter lollipop plots:\n",
    "  Plots/Regional_parameters.png\n",
    "  Plots/Regional_parameters.pdf\n")


# ==============================================================================
# SECTION 14: SAVE OBJECTS FOR BOOTSTRAP
# ==============================================================================

save(
  final_results,
  regions_present,
  island_width_k,
  K_VALUES,
  file = "Bootstrap_Required_Objects.RData"
)

cat("\n============================================================\n")
cat("BOOTSTRAP INPUT WORKSPACE SAVED\n")
cat("============================================================\n")
cat("File: Bootstrap_Required_Objects.RData\n")
cat("Contains:\n")
cat("  final_results\n")
cat("  regions_present\n")
cat("  island_width_k\n")
cat("  K_VALUES\n")
cat("============================================================\n")


