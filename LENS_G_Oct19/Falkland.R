# ==============================================================================
# MODEL D RECHARGE RATES FOR FALKLAND COMPARISON
# ==============================================================================

falkland_comparison <- final_results %>%
  mutate(
    recharge_fraction_D = rf_D[paste0("rf_", Region)],
    recharge_D_m = Ensemble.Precip * recharge_fraction_D,
    recharge_D_mm = recharge_D_m * 1000
  ) %>%
  select(
    Region,
    Island,
    Ensemble.Precip,
    recharge_fraction_D,
    recharge_D_mm
  )

print(falkland_comparison)

falkland_regional <- final_results %>%
  mutate(
    recharge_fraction_D = rf_D[paste0("rf_", Region)],
    recharge_D_mm = Ensemble.Precip *
      recharge_fraction_D * 1000
  ) %>%
  group_by(Region) %>%
  summarise(
    n = sum(is.finite(recharge_D_mm)),
    mean_rainfall_mm = mean(Ensemble.Precip * 1000, na.rm = TRUE),
    median_rainfall_mm = median(Ensemble.Precip * 1000, na.rm = TRUE),
    recharge_fraction_D = first(recharge_fraction_D),
    mean_recharge_D_mm = mean(recharge_D_mm, na.rm = TRUE),
    median_recharge_D_mm = median(recharge_D_mm, na.rm = TRUE)
  )

print(falkland_regional)


# ==============================================================================
# CHAPMAN–FALKLAND RAINFALL–RECHARGE RELATIONSHIP
# Approximate coordinates digitized from Falkland (1990), Fig. 7
# ==============================================================================

falkland_curve <- data.frame(
  rainfall_mm = c(
    700, 800, 900, 1000, 1100, 1200, 1300, 1400,
    1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300
  ),
  recharge_mm = c(
    140, 170, 200, 230, 270, 310, 350, 395,
    435, 480, 530, 590, 660, 750, 850, 980, 1170
  )
)

falkland_comparison <- final_results %>%
  mutate(
    recharge_fraction_D = rf_D[paste0("rf_", Region)],
    recharge_D_m = Ensemble.Precip * recharge_fraction_D,
    recharge_D_mm = recharge_D_m * 1000,
    
    rainfall_mm = Ensemble.Precip * 1000,
    
    falkland_expected_recharge_mm = approx(
      x = falkland_curve$rainfall_mm,
      y = falkland_curve$recharge_mm,
      xout = rainfall_mm,
      rule = 1
    )$y,
    
    recharge_excess_mm = recharge_D_mm - falkland_expected_recharge_mm,
    
    recharge_ratio = recharge_D_mm / falkland_expected_recharge_mm
  ) %>%
  select(
    Region,
    Island,
    rainfall_mm,
    recharge_fraction_D,
    recharge_D_mm,
    falkland_expected_recharge_mm,
    recharge_excess_mm,
    recharge_ratio
  )

falkland_regional_comparison <- falkland_comparison %>%
  group_by(Region) %>%
  summarise(
    n = sum(!is.na(falkland_expected_recharge_mm)),
    mean_rainfall_mm = mean(rainfall_mm, na.rm = TRUE),
    mean_recharge_D_mm = mean(recharge_D_mm, na.rm = TRUE),
    mean_falkland_recharge_mm = mean(
      falkland_expected_recharge_mm,
      na.rm = TRUE
    ),
    mean_recharge_ratio = mean(
      recharge_ratio,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(falkland_regional_comparison)

