#### Runs the algebraic Framework on the 47 atoll islands
#install.packages("lme4")

install.packages("lme4")
install.packages("DescTools")
install.packages("patchwork")

library(dplyr)
library(ggplot2)
library(purrr)
library(plyr)
library(lme4)
library(DescTools)
library(patchwork)

setwd("C:/Users/Owner/OneDrive/Desktop/Github")  # off-site path directory


# CREATING CURVES FOR Y, B, AND C PARAMETERS
#===============================================================================

# points for b factor (based on figures from Bailey's work)
df_b <- data.frame(
  island_width = c(
    0, 25, 70, 120, 190, 250, 290, 400, 500, 700,
    800, 900, 950, 1000, 1050, 1100, 1150, 1200,
    1500, 1800, 2000, 2200, 2300, 2500, 2600
  ),
  b_factor = c(
    6.9, 5.2, 4, 3, 2, 1.3, 1, 0.6, 0.3, 0.25,
    0.15, 0.1, 0.08, 0.06, 0.045, 0.032, 0.021, 0.015,
    0.010, 0.006, 0.004, 0.003, 0.001, 0.0001, 0.0001
  )
)


# points for c factor
df_c <- data.frame(
  island_width = c(
    0, 125, 250, 275, 500, 750, 800, 900, 1000,
    1050, 1200, 1500, 1600, 2000, 2200, 2500, 2600
  ),
  c_factor = c(
    0.853, 0.885, 0.91, 0.92, 0.955, 0.98, 0.985,
    0.99, 0.995, 0.996, 1, 1, 1, 1, 1, 1, 1
  )
)


# points for y factor
df_y <- data.frame(
  y_factor_widths = c(0, 500, 1000, 1500, 2000, 2500),
  y_factor = c(-1, -9.9, -14.9, -17.8, -19.9, -21.2)
)


# points for k factor
island_width_k <- data.frame(
  width = c(
    150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800,
    850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1300, 1400, 1500, 1600,
    1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500
  ),
  k25 = c(
    1.6, 1.6, 1.58, 1.54, 1.47, 1.36, 1.24, 1.13, 1.07, 1.03, 1.02, 1.01,
    1.01, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
  ),
  k50 = c(
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
  ),
  k75 = c(
    0.74, 0.75, 0.77, 0.79, 0.81, 0.83, 0.86, 0.88, 0.91, 0.93, 0.94, 0.96,
    0.97, 0.98, 0.98, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99,
    0.99, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
  ),
  k100 = c(
    0.61, 0.62, 0.63, 0.66, 0.68, 0.72, 0.75, 0.79, 0.83, 0.87, 0.9, 0.92,
    0.94, 0.95, 0.96, 0.97, 0.98, 0.98, 0.98, 0.98, 0.99, 0.99, 0.99, 0.99,
    0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99
  ),
  k125 = c(
    0.52, 0.53, 0.55, 0.57, 0.6, 0.63, 0.67, 0.71, 0.76, 0.8, 0.84, 0.87,
    0.9, 0.92, 0.94, 0.95, 0.96, 0.97, 0.97, 0.97, 0.98, 0.98, 0.98, 0.98,
    0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98
  ),
  k150 = c(
    0.46, 0.47, 0.49, 0.51, 0.53, 0.57, 0.6, 0.65, 0.69, 0.74, 0.78, 0.82,
    0.86, 0.88, 0.91, 0.93, 0.94, 0.95, 0.96, 0.96, 0.97, 0.97, 0.97, 0.98,
    0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98
  ),
  k175 = c(
    0.41, 0.42, 0.44, 0.46, 0.48, 0.51, 0.55, 0.59, 0.63, 0.68, 0.72, 0.77,
    0.81, 0.84, 0.87, 0.89, 0.91, 0.93, 0.94, 0.95, 0.96, 0.96, 0.97, 0.97,
    0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97
  ),
  k200 = c(
    0.38, 0.39, 0.4, 0.42, 0.44, 0.47, 0.5, 0.53, 0.57, 0.62, 0.66, 0.71,
    0.75, 0.79, 0.83, 0.86, 0.88, 0.9, 0.92, 0.93, 0.94, 0.95, 0.96, 0.96,
    0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97
  ),
  k225 = c(
    0.35, 0.36, 0.37, 0.39, 0.4, 0.43, 0.45, 0.49, 0.53, 0.57, 0.61, 0.66,
    0.7, 0.74, 0.78, 0.82, 0.85, 0.87, 0.89, 0.91, 0.92, 0.93, 0.94, 0.95,
    0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96
  ),
  k250 = c(
    0.33, 0.33, 0.34, 0.36, 0.37, 0.4, 0.42, 0.45, 0.48, 0.52, 0.56, 0.61,
    0.65, 0.7, 0.74, 0.77, 0.81, 0.84, 0.86, 0.88, 0.9, 0.91, 0.93, 0.94,
    0.95, 0.95, 0.95, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96
  ),
  k275 = c(
    0.31, 0.31, 0.32, 0.34, 0.35, 0.37, 0.39, 0.42, 0.45, 0.48, 0.52, 0.56,
    0.61, 0.65, 0.69, 0.73, 0.77, 0.8, 0.83, 0.85, 0.87, 0.89, 0.91, 0.93,
    0.94, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.96, 0.96, 0.96
  ),
  k300 = c(
    0.29, 0.3, 0.31, 0.32, 0.33, 0.35, 0.37, 0.39, 0.42, 0.45, 0.49, 0.52,
    0.56, 0.61, 0.65, 0.69, 0.73, 0.76, 0.79, 0.82, 0.85, 0.87, 0.9, 0.92,
    0.94, 0.94, 0.95, 0.95, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96
  ),
  k325 = c(
    0.27, 0.28, 0.29, 0.3, 0.31, 0.33, 0.35, 0.37, 0.4, 0.42, 0.46, 0.49,
    0.53, 0.57, 0.61, 0.65, 0.69, 0.73, 0.76, 0.79, 0.82, 0.85, 0.89, 0.92,
    0.94, 0.95, 0.96, 0.97, 0.97, 0.97, 0.98, 0.98, 0.98, 0.98, 0.98
  ),
  k400 = c(
    0.24, 0.25, 0.25, 0.26, 0.27, 0.29, 0.3, 0.32, 0.34, 0.37, 0.39, 0.42,
    0.45, 0.49, 0.53, 0.56, 0.6, 0.64, 0.68, 0.71, 0.74, 0.77, 0.82, 0.86,
    0.89, 0.91, 0.92, 0.93, 0.94, 0.94, 0.95, 0.95, 0.95, 0.95, 0.95
  )
)


ggplot(island_width_k, aes(x = width, y = k25)) +
  geom_line() +
  geom_line(aes(x = width, y = k50)) +
  geom_line(aes(x = width, y = k75)) +
  geom_line(aes(x = width, y = k100)) +
  geom_line(aes(x = width, y = k125)) +
  geom_line(aes(x = width, y = k150)) +
  geom_line(aes(x = width, y = k175)) +
  geom_line(aes(x = width, y = k200)) +
  geom_line(aes(x = width, y = k225)) +
  geom_line(aes(x = width, y = k250)) +
  geom_line(aes(x = width, y = k275)) +
  geom_line(aes(x = width, y = k300)) +
  geom_line(aes(x = width, y = k325)) +
  geom_line(aes(x = width, y = k400))


write.csv(island_width_k, "Data_Files/island_width_k.csv")


#-------------------------------------------------------------------------------
# DEFINING HELPER FUNCTIONS
#-------------------------------------------------------------------------------

# NOTE:
# The precipitation-product *_R columns in FOR_PROJECTION.csv
# are already recharge values in meters/year.
#
# Therefore, no additional conversion from precipitation to recharge
# is performed in this script.

# C factor
get_C <- function(width) {
  1 - 0.15 * exp(-0.00288 * width)
}


curve_c <- data.frame(island_width = seq(0, 2600, by = 10))
curve_c$fit <- get_C(curve_c$island_width)

ggplot() +
  geom_point(data = df_c, aes(island_width, c_factor), size = 2) +
  geom_line(
    data = curve_c,
    aes(island_width, fit),
    color = "red",
    linewidth = 1.2
  ) +
  labs(
    x = "Island Width (m)",
    y = "C Factor",
    title = "C Scaling Curve"
  ) +
  theme_minimal()


# B factor
get_B <- function(width) {
  5 * exp(-0.00461 * width)
}


curve_b <- data.frame(island_width = seq(0, 2600, by = 10))
curve_b$fit <- get_B(curve_b$island_width)

ggplot() +
  geom_point(data = df_b, aes(island_width, b_factor), size = 2) +
  geom_line(
    data = curve_b,
    aes(island_width, fit),
    color = "red",
    linewidth = 1.2
  ) +
  labs(
    x = "Island Width (m)",
    y = "B Factor",
    title = "B Scaling Curve"
  ) +
  theme_minimal()


# Y factor
get_Y <- function(width) {
  -20 * (1 - exp(-0.001386 * width))
}


curve_y <- data.frame(width = seq(0, 2600, by = 10))
curve_y$fit <- get_Y(curve_y$width)

ggplot() +
  geom_point(data = df_y, aes(y_factor_widths, y_factor), size = 2) +
  geom_line(
    data = curve_y,
    aes(width, fit),
    color = "red",
    linewidth = 1.2
  ) +
  labs(
    x = "Island Width (m)",
    y = "Y Factor",
    title = "Y Scaling Curve"
  ) +
  theme_minimal()


# Return K parameter; rounding differently based on how wide the island is
get_K <- function(width, hc) {
  
  # Handle missing or invalid widths
  if (is.na(width) || !is.numeric(width)) {
    print("Width is NA or non-numeric; returning NA for K")
    return(NA_real_)
  }
  
  # Rounding logic
  if (width < 101) {
    width_rounded <- 150
  } else if (width < 1249) {
    width_rounded <- round_any(width, 50)
  } else if (width < 2500) {
    width_rounded <- round_any(width, 100)
  } else {
    width_rounded <- 2500
  }
  
  print(paste("Rounded Width", width_rounded, sep = ": "))
  
  # Validate lookup row
  K_row <- which(island_width_k$width == width_rounded)
  
  if (length(K_row) == 0) {
    print(paste(
      "No matching width in lookup table for",
      width_rounded
    ))
    return(NA_real_)
  }
  
  # Validate lookup column
  K_column <- which(colnames(island_width_k) == hc)
  
  if (length(K_column) == 0) {
    print(paste(
      "No matching K column for",
      hc
    ))
    return(NA_real_)
  }
  
  # Extract K
  K <- island_width_k[K_row, K_column]
  
  return(K)
}


#-------------------------------------------------------------------------------
# DEFINING AND RUNNING THE ALGEBRAIC MODEL
#-------------------------------------------------------------------------------

Bailey_Model <- function(Y, Z_TD, R, B, K, C) {
  
  Z_MAX = ((((Z_TD - Y) * R) / (B + R)) + Y) * K * C
  
  return(as.numeric(Z_MAX))
}


# Running all functions to obtain parameters for each island
get_Z_MAX <- function(R, width, hc, Z_TD) {
  
  B <- get_B(width)
  C <- get_C(width)
  Y <- get_Y(width)
  K <- get_K(width, hc)
  
  Z_MAX <- Bailey_Model(
    Y,
    Z_TD,
    R,
    B,
    K,
    C
  )
  
  return(as.numeric(Z_MAX))
}


#-------------------------------------------------------------------------------
# READ INPUT DATA AND CALCULATE ENSEMBLE RECHARGE
#-------------------------------------------------------------------------------

island_features <- read.csv("Data_Files/FOR_PROJECTION.csv")


# The five precipitation products have already been converted to recharge
# in meters/year:
#
#   MSWEP_R
#   JAXA_R
#   NASA_R
#   CDR_R
#   ERA5_R
#
# Calculate the five-product ensemble recharge for every island.
#
# na.rm = FALSE means that all five products must be available for an island
# to receive an ensemble value.

island_features <- island_features %>%
  mutate(
    Ensemble.R = rowMeans(
      cbind(
        MSWEP_R,
        JAXA_R,
        NASA_R,
        CDR_R,
        ERA5_R
      ),
      na.rm = FALSE
    )
  )


# Optional check of the new ensemble recharge values
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
      Ensemble.R
    ) %>%
    head()
)


#-------------------------------------------------------------------------------
# RUN BAILEY MODEL USING ENSEMBLE RECHARGE
#-------------------------------------------------------------------------------

final_results <- island_features %>%
  mutate(
    
    # Central/high-confidence estimate
    Z_MAX_estimated_hc = pmap_dbl(
      list(
        Ensemble.R,
        Avg.Width,
        Estimated.K,
        Z_TD
      ),
      ~ get_Z_MAX(..1, ..2, ..3, ..4)
    ),
    
    # Low estimate
    Z_MAX_estimated_LOW = pmap_dbl(
      list(
        Ensemble.R,
        Width.1,
        High.K,
        Z_TD
      ),
      ~ get_Z_MAX(..1, ..2, ..3, ..4)
    ),
    
    # High estimate
    Z_MAX_estimated_HIGH = pmap_dbl(
      list(
        Ensemble.R,
        Width.4,
        Low.K,
        Z_TD
      ),
      ~ get_Z_MAX(..1, ..2, ..3, ..4)
    )
  )


# Save results
write.csv(
  final_results,
  "Data_Files/Bailey_results.csv"
)


#-------------------------------------------------------------------------------
# MARSHALL ISLANDS AND FSM
#-------------------------------------------------------------------------------

final_results2 <- final_results %>%
  filter(
    Region == "FSM" |
      Region == "Marshalls"
  ) %>%
  mutate(
    residual = Z_MAX_estimated_hc - Z_MAX_O
  )


# WITHOUT BAILEY'S PREDICTIONS
ggplot(
  final_results2,
  aes(
    x = Z_MAX_O,
    y = Z_MAX_estimated_hc,
    shape = Region
  )
) +
  geom_abline(
    aes(color = "1:1 line"),
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 1,
    color = "red"
  ) +
  geom_point(size = 3) +
  xlim(0, 25) +
  ylim(0, 25) +
  labs(
    x = "Observed Maximum Lens Thickness (m)",
    y = "Predicted Maximum Lens Thickness (m)"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "right",
    legend.text = element_text(size = 12)
  )


ggsave(
  filename = "R_Plots/MARSH_FSM.png",
  width = 7.5,
  height = 6,
  dpi = 300
)


# WITH BAILEY'S PREDICTIONS
ggplot(
  final_results2,
  aes(
    x = Z_MAX_O,
    y = Z_MAX_estimated_hc,
    shape = Region
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 1,
    color = "red"
  ) +
  geom_point(
    size = 3,
    aes(color = "Our Model")
  ) +
  geom_point(
    aes(
      y = Z_MAX_B,
      color = "Bailey's Predictions"
    ),
    size = 3
  ) +
  lims(
    x = c(0, 25),
    y = c(0, 25)
  ) +
  scale_color_manual(
    name = "",
    values = c(
      "1:1 line" = "red",
      "Our Model" = "black",
      "Bailey's Predictions" = "cyan"
    )
  ) +
  labs(
    x = "Observed Maximum Lens Thickness (m)",
    y = "Predicted Maximum Lens Thickness (m)"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "right",
    legend.text = element_text(size = 12)
  )


ggsave(
  filename = "R_Plots/MARSH_FSM_BAILEY.png",
  width = 8,
  height = 6,
  dpi = 300
)


#-------------------------------------------------------------------------------
# PLOTTING RESIDUALS
#-------------------------------------------------------------------------------

ggplot(
  final_results2,
  aes(
    x = Z_MAX_O,
    y = residual
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "red"
  ) +
  geom_point(size = 2.5) +
  lims(y = c(-10, 10)) +
  labs(
    x = "Observed Maximum Lens Thickness (m)",
    y = "Residual"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )


ggsave(
  filename = "R_Plots/MARSH_FSM_RES.png",
  width = 7,
  height = 8,
  dpi = 300
)


#-------------------------------------------------------------------------------
# QUANTIFYING THE FIT OF THE DATA
#-------------------------------------------------------------------------------

# RMSE
RMSE <- sqrt(
  mean(
    (final_results2$Z_MAX_estimated_hc -
       final_results2$Z_MAX_O)^2,
    na.rm = TRUE
  )
)

RMSE


# BIAS
bias <- mean(
  final_results2$Z_MAX_estimated_hc -
    final_results2$Z_MAX_O,
  na.rm = TRUE
)

bias


# MAE
MAE <- mean(
  abs(
    final_results2$Z_MAX_estimated_hc -
      final_results2$Z_MAX_O
  ),
  na.rm = TRUE
)

MAE


# CCC
CCC <- CCC(
  final_results2$Z_MAX_O,
  final_results2$Z_MAX_estimated_hc
)$rho.c

CCC


# R2 from regressing predicted on observed with slope forced to 1
model_force <- lm(
  Z_MAX_estimated_hc - Z_MAX_O ~ 0,
  data = final_results2
)

summary(model_force)


#-------------------------------------------------------------------------------
# ALL REGIONS
#-------------------------------------------------------------------------------

final_results <- final_results %>%
  filter(
    Region == "Marshalls" |
      Region == "FSM" |
      Region == "Maldives" |
      Region == "Central_Indian" |
      Region == "East_Indian" |
      Region == "Kiribati" |
      Region == "North_Pacific" |
      Region == "Cook_Islands"
  )


# Define exact region order
final_results$Region <- factor(
  final_results$Region,
  levels = c(
    "Marshalls",
    "FSM",
    "Maldives",
    "Gilbert",
    "Central_Indian",
    "East_Indian",
    "Kiribati",
    "North_Pacific",
    "Cook_Islands",
    "Indonesia"
  )
)


# Define colors
region_colors <- c(
  "Marshalls" = "firebrick",
  "FSM" = "cyan2",
  "Maldives" = "darkorchid",
  "Gilbert" = "green",
  "Central_Indian" = "magenta",
  "East_Indian" = "orange",
  "Kiribati" = "darkcyan",
  "North_Pacific" = "blue",
  "Cook_Islands" = "maroon"
)


ggplot(
  final_results,
  aes(
    x = Z_MAX_O,
    y = Z_MAX_estimated_hc,
    col = Region
  )
) +
  geom_abline(
    aes(color = "1:1 line"),
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 1,
    color = "red"
  ) +
  geom_point(size = 3.5) +
  xlim(0, 30) +
  ylim(0, 30) +
  scale_color_manual(
    values = region_colors,
    labels = c(
      "Marshall Islands",
      "Micronesia",
      "Maldives",
      "Central Indian",
      "East Indian",
      "Kiribati",
      "North Pacific",
      "Cook Islands"
    )
  ) +
  labs(
    x = "Observed Lens Thickness (m)",
    y = "Predicted Lens Thickness (m)"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 11)
  )


ggsave(
  filename = "R_Plots/All_Regions.png",
  width = 9,
  height = 7,
  dpi = 300
)

