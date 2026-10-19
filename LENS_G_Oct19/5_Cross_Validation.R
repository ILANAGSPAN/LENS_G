load("C:/Users/Owner/OneDrive/Desktop/Github/Bootstrap_Required_Objects.RData")

# ==============================================================================
# Cross_Validation_fast.R
# Transfer validation: 4 largest regions -> 4 smallest regions
# Fast/vectorized implementation with 10,000-replicate bootstrap
# ==============================================================================

setwd("C:/Users/Owner/OneDrive/Desktop/Github/R_Scripts")

library(parallel)
library(DEoptim)
library(dplyr)

# ------------------------------------------------------------------------------
# 1. Required objects
# ------------------------------------------------------------------------------

required_objects <- c(
  "final_results", "regions_present", "n_regions", "island_width_k",
  "K_VALUES", "get_K_orientation", "get_K_correction", "bailey_core"
)

missing <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing) > 0) {
  stop(paste(
    "Missing required objects:", paste(missing, collapse = ", "),
    "\nRun Optimization_Test.R first."
  ))
}

# ------------------------------------------------------------------------------
# 2. Dataset and calibration/validation split
# ------------------------------------------------------------------------------

model_data <- final_results %>%
  filter(
    is.finite(Ensemble.Precip),
    is.finite(Avg.Width),
    is.finite(Z_TD),
    is.finite(Z_MAX_O),
    !is.na(Orientation),
    !is.na(Region)
  )

cat("\n============================================================\n")
cat("TRANSFER VALIDATION WITH BOOTSTRAP UNCERTAINTY\n")
cat("FAST / VECTORIZED IMPLEMENTATION\n")
cat("============================================================\n\n")

cat("Total usable observations:", nrow(model_data), "\n")
cat("Regions:", n_distinct(model_data$Region), "\n\n")

region_counts <- model_data %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n))

print(region_counts)

cat("Observations by region:\n\n")
print(region_counts)

large_regions <- region_counts %>%
  slice_head(n = 4) %>%
  pull(Region)

small_regions <- region_counts %>%
  slice_tail(n = 4) %>%
  pull(Region)

cat("\nFour largest regions (CALIBRATION):\n")
print(region_counts %>% filter(Region %in% large_regions))

cat("\nFour smallest regions (VALIDATION):\n")
print(region_counts %>% filter(Region %in% small_regions))

large_data <- model_data %>%
  filter(Region %in% large_regions)

small_data <- model_data %>%
  filter(Region %in% small_regions)

cat("\nCalibration observations:", nrow(large_data), "\n")
cat("Validation observations: ", nrow(small_data), "\n\n")

# Stable ID for precomputed K corrections
large_data$.transfer_id <- seq_len(nrow(large_data))

# ------------------------------------------------------------------------------
# 3. Precompute K-correction lookup
# ------------------------------------------------------------------------------

width_clamped_large <- pmin(pmax(large_data$Avg.Width, 150), 2500)

Kcorr_large <- sapply(
  K_VALUES,
  function(K_value) {
    approx(
      island_width_k$width,
      island_width_k[[paste0("k", K_value)]],
      xout = width_clamped_large,
      rule = 2
    )$y
  }
)

Kcorr_large <- as.matrix(Kcorr_large)
colnames(Kcorr_large) <- paste0(".Kcorr_", K_VALUES)

cat("Precomputed K-correction matrix:\n")
cat("  Calibration islands:", nrow(Kcorr_large), "\n")
cat("  K values:", ncol(Kcorr_large), "\n\n")

# ------------------------------------------------------------------------------
# 4. Original slow objective, retained for verification only
# ------------------------------------------------------------------------------

obj_fn_transfer_slow <- function(params, data) {
  rf_global <- params[1]
  Kscale_global <- params[2]
  
  pred <- mapply(
    FUN = function(rain, width, ztd, orient) {
      K_mday <- get_K_orientation(orient) * Kscale_global
      K_mday <- pmin(pmax(K_mday, 25), 400)
      
      bailey_core(
        rainfall_m = rain,
        width = width,
        Z_TD = ztd,
        recharge_frac = rf_global,
        K_mday = K_mday
      )
    },
    data$Ensemble.Precip,
    data$Avg.Width,
    data$Z_TD,
    data$Orientation
  )
  
  valid <- is.finite(pred) & is.finite(data$Z_MAX_O)
  
  if (sum(valid) < 10) return(1e10)
  
  sqrt(mean((pred[valid] - data$Z_MAX_O[valid])^2))
}

# ------------------------------------------------------------------------------
# 5. Fast vectorized objective
# ------------------------------------------------------------------------------

obj_fn_transfer_fast <- function(params, data, Kcorr_matrix, K_values) {
  rf_global <- params[1]
  Kscale_global <- params[2]
  
  orientation <- as.character(data$Orientation)
  
  K_orientation <- rep(200, nrow(data))
  K_orientation[!is.na(orientation) & orientation == "windward"] <- 400
  K_orientation[!is.na(orientation) & orientation == "leeward"] <- 50
  K_orientation[!is.na(orientation) & orientation == "offset"] <- 200
  
  K_mday <- K_orientation * Kscale_global
  K_mday <- pmin(pmax(K_mday, 25), 400)
  K_mday <- pmin(pmax(K_mday, min(K_values)), max(K_values))
  
  lower_index <- findInterval(K_mday, K_values)
  lower_index <- pmin(pmax(lower_index, 1), length(K_values) - 1)
  
  upper_index <- lower_index + 1
  
  upper_index[K_mday >= max(K_values)] <- length(K_values)
  lower_index[K_mday >= max(K_values)] <- length(K_values)
  
  row_index <- data$.transfer_id
  
  correction_lower <- Kcorr_matrix[
    cbind(row_index, lower_index)
  ]
  
  correction_upper <- Kcorr_matrix[
    cbind(row_index, upper_index)
  ]
  
  K_lower <- K_values[lower_index]
  K_upper <- K_values[upper_index]
  
  K_correction <- ifelse(
    K_upper == K_lower,
    correction_lower,
    correction_lower +
      ((K_mday - K_lower) / (K_upper - K_lower)) *
      (correction_upper - correction_lower)
  )
  
  rainfall_m <- data$Ensemble.Precip
  width <- data$Avg.Width
  Z_TD <- data$Z_TD
  
  R <- rainfall_m * rf_global
  B <- 5 * exp(-0.00461 * width)
  C <- 1 - 0.15 * exp(-0.00288 * width)
  Y <- -20 * (1 - exp(-0.001386 * width))
  
  Z_MAX <- (
    ((Z_TD - Y) * R) / (B + R) + Y
  ) * K_correction * C
  
  valid <- is.finite(Z_MAX) & is.finite(data$Z_MAX_O)
  
  if (sum(valid) < 10) return(1e10)
  
  sqrt(mean((Z_MAX[valid] - data$Z_MAX_O[valid])^2))
}

# ------------------------------------------------------------------------------
# 6. Verify fast objective against original objective
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("VERIFYING FAST OBJECTIVE AGAINST ORIGINAL OBJECTIVE\n")
cat("============================================================\n\n")

set.seed(12345)

n_objective_tests <- 25

test_params <- data.frame(
  rf_global = runif(n_objective_tests, 0.10, 0.90),
  Kscale_global = runif(n_objective_tests, 0.25, 4.00)
)

slow_values <- numeric(n_objective_tests)
fast_values <- numeric(n_objective_tests)

for (i in seq_len(n_objective_tests)) {
  params_i <- c(
    test_params$rf_global[i],
    test_params$Kscale_global[i]
  )
  
  slow_values[i] <- obj_fn_transfer_slow(
    params_i,
    large_data
  )
  
  fast_values[i] <- obj_fn_transfer_fast(
    params_i,
    large_data,
    Kcorr_large,
    K_VALUES
  )
}

objective_differences <- abs(slow_values - fast_values)
max_objective_difference <- max(objective_differences, na.rm = TRUE)
objective_tolerance <- 1e-8

cat(
  "Number of parameter combinations tested:",
  n_objective_tests,
  "\n"
)

cat(
  "Maximum absolute difference:",
  format(max_objective_difference, scientific = TRUE),
  "\n\n"
)

if (
  !is.finite(max_objective_difference) ||
  max_objective_difference > objective_tolerance
) {
  verification_table <- data.frame(
    rf_global = test_params$rf_global,
    Kscale_global = test_params$Kscale_global,
    slow = slow_values,
    fast = fast_values,
    absolute_difference = objective_differences
  )
  
  print(verification_table)
  
  stop(paste0(
    "\nThe fast transfer objective does not match the original ",
    "objective within the required tolerance of ",
    objective_tolerance,
    ".\nThe bootstrap has NOT been started."
  ))
}

cat("OBJECTIVE VERIFICATION PASSED.\n")
cat(
  "The fast and original objectives agree within ",
  objective_tolerance,
  ".\n\n"
)

objective_verification <- data.frame(
  rf_global = test_params$rf_global,
  Kscale_global = test_params$Kscale_global,
  slow = slow_values,
  fast = fast_values,
  absolute_difference = objective_differences
)

write.csv(
  objective_verification,
  "transfer_objective_verification.csv",
  row.names = FALSE
)

cat("Saved: transfer_objective_verification.csv\n\n")

# ------------------------------------------------------------------------------
# 7. Bootstrap settings
# ------------------------------------------------------------------------------

n_boot <- 10000

available_cores <- parallel::detectCores(logical = TRUE)

if (is.na(available_cores)) {
  n_cores <- 1
} else {
  n_cores <- max(1, available_cores - 1)
}

bootstrap_itermax <- 50
bootstrap_steptol <- 10
bootstrap_reltol <- 1e-5
bootstrap_CR <- 0.9
bootstrap_F <- 0.8
batch_size <- 500
use_initial_population <- FALSE
bootstrap_NP <- 8

cat("============================================================\n")
cat("BOOTSTRAP SETTINGS\n")
cat("============================================================\n")
cat("Bootstrap replicates:        ", n_boot, "\n")
cat("Logical CPU cores available: ", available_cores, "\n")
cat("Cores used:                  ", n_cores, "\n")
cat("DEoptim NP:                  ", bootstrap_NP, "\n")
cat("DEoptim itermax:             ", bootstrap_itermax, "\n")
cat("DEoptim steptol:             ", bootstrap_steptol, "\n")
cat("DEoptim reltol:              ", bootstrap_reltol, "\n")
cat("DEoptim CR:                  ", bootstrap_CR, "\n")
cat("DEoptim F:                   ", bootstrap_F, "\n")
cat("Batch size:                  ", batch_size, "\n")
cat("============================================================\n\n")

# ------------------------------------------------------------------------------
# 8. One bootstrap replicate
# ------------------------------------------------------------------------------

one_bootstrap_replicate_transfer <- function(
    rep_idx,
    calib_data,
    boot_itermax,
    boot_NP,
    boot_steptol,
    boot_reltol,
    boot_CR,
    boot_F,
    use_init_pop,
    Kcorr_matrix,
    K_values
) {
  idx <- sample(
    seq_len(nrow(calib_data)),
    size = nrow(calib_data),
    replace = TRUE
  )
  
  boot_data <- calib_data[idx, , drop = FALSE]
  rownames(boot_data) <- NULL
  
  lower <- c(0.10, 0.25)
  upper <- c(0.90, 4.00)
  
  init_pop <- NULL
  
  if (use_init_pop) {
    init_pop <- matrix(
      runif(boot_NP * 2, min = lower, max = upper),
      nrow = boot_NP,
      ncol = 2
    )
  }
  
  fit <- tryCatch(
    DEoptim(
      fn = obj_fn_transfer_fast,
      lower = lower,
      upper = upper,
      data = boot_data,
      Kcorr_matrix = Kcorr_matrix,
      K_values = K_values,
      control = DEoptim.control(
        itermax = boot_itermax,
        NP = boot_NP,
        CR = boot_CR,
        F = boot_F,
        trace = TRUE,
        steptol = boot_steptol,
        reltol = boot_reltol,
        storepopfrom = boot_itermax + 1,
        initialpop = init_pop
      )
    ),
    error = function(e) {
      list(.error = conditionMessage(e))
    }
  )
  
  if (is.list(fit) && ".error" %in% names(fit)) {
    return(list(
      success = FALSE,
      replicate = rep_idx,
      params = NULL,
      objective = NA_real_,
      error = fit$.error
    ))
  }
  
  result <- c(
    rf_global = fit$optim$bestmem[1],
    Kscale_global = fit$optim$bestmem[2]
  )
  
  list(
    success = TRUE,
    replicate = rep_idx,
    params = result,
    objective = fit$optim$bestval,
    error = NA_character_
  )
}

# ------------------------------------------------------------------------------
# 9. Initialize parallel cluster
# ------------------------------------------------------------------------------

cat("\n")
cat("INITIALIZING PARALLEL BOOTSTRAP\n\n")

cat(
  "Creating PSOCK cluster with ",
  n_cores,
  " workers...\n",
  sep = ""
)

cl <- parallel::makeCluster(
  n_cores,
  type = "PSOCK"
)

parallel::clusterEvalQ(
  cl,
  {
    library(DEoptim)
    library(dplyr)
    NULL
  }
)

objects_to_export <- c(
  "one_bootstrap_replicate_transfer",
  "obj_fn_transfer_fast",
  "large_data",
  "Kcorr_large",
  "K_VALUES",
  "bootstrap_itermax",
  "bootstrap_NP",
  "bootstrap_steptol",
  "bootstrap_reltol",
  "bootstrap_CR",
  "bootstrap_F",
  "use_initial_population"
)

parallel::clusterExport(
  cl,
  varlist = objects_to_export,
  envir = .GlobalEnv
)

parallel::clusterSetRNGStream(cl)

worker_test <- parallel::clusterEvalQ(
  cl,
  {
    list(
      has_bootstrap_function =
        exists("one_bootstrap_replicate_transfer"),
      has_fast_objective =
        exists("obj_fn_transfer_fast"),
      has_Kcorr_matrix =
        exists("Kcorr_large"),
      has_K_values =
        exists("K_VALUES")
    )
  }
)

print(worker_test)

worker_test_ok <- all(
  vapply(
    worker_test,
    function(x) all(unlist(x)),
    logical(1)
  )
)

if (!worker_test_ok) {
  parallel::stopCluster(cl)
  stop("One or more workers failed the initialization test.")
}

cat("\nAll workers passed initialization test.\n\n")

# ------------------------------------------------------------------------------
# 10. Test one bootstrap replicate
# ------------------------------------------------------------------------------

cat("Testing one bootstrap replicate before full run...\n")

test_result <- parallel::parLapply(
  cl,
  1,
  function(i) {
    one_bootstrap_replicate_transfer(
      i,
      large_data,
      bootstrap_itermax,
      bootstrap_NP,
      bootstrap_steptol,
      bootstrap_reltol,
      bootstrap_CR,
      bootstrap_F,
      use_initial_population,
      Kcorr_large,
      K_VALUES
    )
  }
)

if (
  length(test_result) != 1 ||
  !isTRUE(test_result[[1]]$success)
) {
  parallel::stopCluster(cl)
  
  stop(
    paste0(
      "The single bootstrap test failed: ",
      test_result[[1]]$error
    )
  )
}

cat("Single bootstrap test passed.\n")
cat(
  "Test RF: ",
  test_result[[1]]$params["rf_global"],
  "\n",
  sep = ""
)
cat(
  "Test K scalar: ",
  test_result[[1]]$params["Kscale_global"],
  "\n\n",
  sep = ""
)

# ------------------------------------------------------------------------------
# 11. Run 10,000 bootstrap replicates
# ------------------------------------------------------------------------------

cat("################################################################################\n")
cat("RUNNING 10,000 BOOTSTRAP REPLICATES\n")
cat("################################################################################\n\n")

start_time <- Sys.time()
boot_results <- vector("list", n_boot)

batch_starts <- seq(
  1,
  n_boot,
  by = batch_size
)

for (batch_start in batch_starts) {
  batch_end <- min(
    batch_start + batch_size - 1,
    n_boot
  )
  
  batch_indices <- batch_start:batch_end
  
  cat(
    "Bootstrap replicates ",
    batch_start,
    " - ",
    batch_end,
    "...\n",
    sep = ""
  )
  
  batch_results <- tryCatch(
    parallel::parLapply(
      cl,
      batch_indices,
      function(i) {
        one_bootstrap_replicate_transfer(
          i,
          large_data,
          bootstrap_itermax,
          bootstrap_NP,
          bootstrap_steptol,
          bootstrap_reltol,
          bootstrap_CR,
          bootstrap_F,
          use_initial_population,
          Kcorr_large,
          K_VALUES
        )
      }
    ),
    error = function(e) {
      stop(
        paste(
          "Parallel bootstrap failed:",
          conditionMessage(e)
        )
      )
    }
  )
  
  boot_results[batch_indices] <- batch_results
  
  completed <- batch_end
  
  successful_so_far <- sum(
    vapply(
      boot_results[seq_len(completed)],
      function(x) {
        !is.null(x) && isTRUE(x$success)
      },
      logical(1)
    )
  )
  
  elapsed_seconds <- as.numeric(
    difftime(
      Sys.time(),
      start_time,
      units = "secs"
    )
  )
  
  rate <- completed / elapsed_seconds
  
  remaining_seconds <- if (rate > 0) {
    (n_boot - completed) / rate
  } else {
    NA_real_
  }
  
  percent_complete <- 100 * completed / n_boot
  
  cat(
    "  Progress: ",
    completed,
    " / ",
    n_boot,
    " (",
    round(percent_complete, 1),
    "%)\n",
    sep = ""
  )
  
  cat(
    "  Successful: ",
    successful_so_far,
    " | Elapsed: ",
    round(elapsed_seconds / 3600, 2),
    " hrs | Est. remaining: ",
    round(remaining_seconds / 3600, 2),
    " hrs\n\n",
    sep = ""
  )
}

parallel::stopCluster(cl)

end_time <- Sys.time()

elapsed_time <- difftime(
  end_time,
  start_time,
  units = "hours"
)

# ------------------------------------------------------------------------------
# 12. Bootstrap diagnostics
# ------------------------------------------------------------------------------

successful <- boot_results[
  vapply(
    boot_results,
    function(x) {
      !is.null(x) && isTRUE(x$success)
    },
    logical(1)
  )
]

n_successful <- length(successful)
success_rate <- 100 * n_successful / n_boot

cat("BOOTSTRAP COMPLETE\n")
cat(
  "Successful: ",
  n_successful,
  " / ",
  n_boot,
  " (",
  round(success_rate, 1),
  "%)\n",
  sep = ""
)
cat(
  "Elapsed time: ",
  round(as.numeric(elapsed_time), 2),
  " hours\n\n",
  sep = ""
)

if (n_successful < n_boot) {
  failed <- boot_results[
    vapply(
      boot_results,
      function(x) {
        !is.null(x) && !isTRUE(x$success)
      },
      logical(1)
    )
  ]
  
  if (length(failed) > 0) {
    cat("\n*** DIAGNOSTIC: Bootstrap failures detected ***\n")
    cat(
      "Number of failures: ",
      length(failed),
      "\n\n",
      sep = ""
    )
    
    cat("First 5 error messages:\n\n")
    
    for (i in seq_len(min(5, length(failed)))) {
      if (!is.na(failed[[i]]$error)) {
        cat(
          "Replicate ",
          failed[[i]]$replicate,
          ": ",
          failed[[i]]$error,
          "\n",
          sep = ""
        )
      }
    }
    
    cat("\n")
  }
}

# ------------------------------------------------------------------------------
# 13. Extract bootstrap parameters
# ------------------------------------------------------------------------------

successful_ids <- vapply(
  successful,
  function(x) x$replicate,
  integer(1)
)

successful_params <- do.call(
  rbind,
  lapply(
    successful,
    function(x) x$params
  )
)

boot_params_with_id <- data.frame(
  replicate = successful_ids,
  rf_global = successful_params[, 1],
  Kscale_global = successful_params[, 2]
)

boot_params_with_id <- boot_params_with_id[
  order(boot_params_with_id$replicate),
]

rownames(boot_params_with_id) <- NULL

boot_params_mat <- boot_params_with_id %>%
  select(rf_global, Kscale_global)

write.csv(
  boot_params_mat,
  "bootstrap_Transfer_parameters.csv",
  row.names = FALSE
)

write.csv(
  boot_params_with_id,
  "bootstrap_Transfer_parameters_with_replicate_id.csv",
  row.names = FALSE
)

cat(
  "Saved: Data_Files/bootstrap_Transfer_parameters.csv (",
  nrow(boot_params_mat),
  " successful replicates)\n",
  sep = ""
)

cat(
  "Saved: Data_Files/bootstrap_Transfer_parameters_with_replicate_id.csv\n\n"
)

# ------------------------------------------------------------------------------
# 14. Select all successful bootstrap replicates
# ------------------------------------------------------------------------------

boot_params_final <- boot_params_with_id

n_final_boot <- nrow(boot_params_final)

if (n_final_boot != n_boot) {
  stop(
    paste0(
      "\nExpected ", n_boot,
      " successful bootstrap replicates, but found ",
      n_final_boot,
      ".\nAll replicates must be successful before calculating ",
      "the 95% bootstrap prediction intervals."
    )
  )
}

cat(
  "Confirmed ", n_final_boot,
  " successful bootstrap replicates for prediction uncertainty.\n\n",
  sep = ""
)

# ------------------------------------------------------------------------------
# 15. Prediction function
# ------------------------------------------------------------------------------

predict_transfer <- function(data, rf_global, Kscale_global) {
  mapply(
    FUN = function(rain, width, ztd, orient) {
      K_mday <- get_K_orientation(orient) * Kscale_global
      K_mday <- pmin(pmax(K_mday, 25), 400)
      
      bailey_core(
        rainfall_m = rain,
        width = width,
        Z_TD = ztd,
        recharge_frac = rf_global,
        K_mday = K_mday
      )
    },
    data$Ensemble.Precip,
    data$Avg.Width,
    data$Z_TD,
    data$Orientation
  )
}

# ------------------------------------------------------------------------------
# 16. Generate validation predictions
# ------------------------------------------------------------------------------

cat(
  "Generating predictions using all 10,000 successful bootstrap replicates...\n\n"
)

n_val <- nrow(small_data)

predictions_matrix <- matrix(
  NA_real_,
  nrow = n_val,
  ncol = nrow(boot_params_final)
)

for (i in seq_len(nrow(boot_params_final))) {
  
  predictions_matrix[, i] <- predict_transfer(
    small_data,
    boot_params_final$rf_global[i],
    boot_params_final$Kscale_global[i]
  )
  
  if (i %% 500 == 0) {
    cat(
      "  Completed predictions for bootstrap replicate ",
      i,
      " / ",
      nrow(boot_params_final),
      "\n",
      sep = ""
    )
  }
}

cat(
  "Predictions complete (10,000 replicates).\n\n"
)

# ------------------------------------------------------------------------------
# 17. 95% Bootstrap Prediction Intervals
# ------------------------------------------------------------------------------

mean_pred <- apply(
  predictions_matrix,
  1,
  mean,
  na.rm = TRUE
)

median_pred <- apply(
  predictions_matrix,
  1,
  median,
  na.rm = TRUE
)

sd_pred <- apply(
  predictions_matrix,
  1,
  sd,
  na.rm = TRUE
)

q025_pred <- apply(
  predictions_matrix,
  1,
  quantile,
  probs = 0.025,
  na.rm = TRUE
)

q975_pred <- apply(
  predictions_matrix,
  1,
  quantile,
  probs = 0.975,
  na.rm = TRUE
)

small_data_with_uncertainty <- small_data %>%
  mutate(
    Z_MAX_mean = mean_pred,
    Z_MAX_median = median_pred,
    Z_MAX_sd = sd_pred,
    Z_MAX_lower_pi_95 = q025_pred,
    Z_MAX_upper_pi_95 = q975_pred,
    Z_MAX_pi_width_95 = q975_pred - q025_pred,
    Z_MAX_error = Z_MAX_mean - Z_MAX_O,
    covered_95 = (
      Z_MAX_O >= Z_MAX_lower_pi_95 &
        Z_MAX_O <= Z_MAX_upper_pi_95
    )
  )

# Overall 95% prediction-interval coverage
coverage_overall <- data.frame(
  n = nrow(small_data_with_uncertainty),
  n_covered = sum(
    small_data_with_uncertainty$covered_95,
    na.rm = TRUE
  ),
  coverage = mean(
    small_data_with_uncertainty$covered_95,
    na.rm = TRUE
  ),
  mean_width = mean(
    small_data_with_uncertainty$Z_MAX_pi_width_95,
    na.rm = TRUE
  ),
  median_width = median(
    small_data_with_uncertainty$Z_MAX_pi_width_95,
    na.rm = TRUE
  )
)

coverage_by_region <- small_data_with_uncertainty %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    n = dplyr::n(),
    n_covered = sum(covered_95, na.rm = TRUE),
    coverage = mean(covered_95, na.rm = TRUE),
    mean_width = mean(
      Z_MAX_pi_width_95,
      na.rm = TRUE
    ),
    median_width = median(
      Z_MAX_pi_width_95,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

cat("============================================================\n")
cat("95% BOOTSTRAP PREDICTION-INTERVAL COVERAGE\n")
cat("============================================================\n\n")

cat("Overall coverage:\n")
print(coverage_overall)

cat("\nCoverage by region:\n")
print(coverage_by_region)

cat("\n")

# ------------------------------------------------------------------------------
# 18. Performance metrics
# ------------------------------------------------------------------------------

calculate_metrics_cv <- function(
    data,
    prediction_column = "Z_MAX_mean"
) {
  pred <- data[[prediction_column]]
  obs <- data$Z_MAX_O
  
  valid <- is.finite(pred) & is.finite(obs)
  
  pred <- pred[valid]
  obs <- obs[valid]
  
  data.frame(
    n = length(obs),
    RMSE = sqrt(mean((pred - obs)^2)),
    MAE = mean(abs(pred - obs)),
    Bias = mean(pred - obs),
    R2 = cor(obs, pred)^2,
    RMSE_std = sd((pred - obs)^2)
  )
}


# ------------------------------------------------------------------------------
# Full-data calibration point estimate
# ------------------------------------------------------------------------------

cat("\nCalibrating transfer model on all calibration observations...\n")

calibration_fit <- DEoptim(
  fn = obj_fn_transfer_fast,
  lower = c(0.10, 0.25),
  upper = c(0.90, 4.00),
  data = large_data,
  Kcorr_matrix = Kcorr_large,
  K_values = K_VALUES,
  control = DEoptim.control(
    itermax = 50,
    NP = 8,
    CR = 0.9,
    F = 0.8,
    trace = FALSE,
    steptol = 10,
    reltol = 1e-5
  )
)

calibration_rf <- calibration_fit$optim$bestmem[1]
calibration_Kscale <- calibration_fit$optim$bestmem[2]

cat(
  "Full-data calibration RF: ",
  calibration_rf,
  "\n",
  sep = ""
)

cat(
  "Full-data calibration K scalar: ",
  calibration_Kscale,
  "\n",
  sep = ""
)

cat(
  "Full-data calibration RMSE from DEoptim: ",
  calibration_fit$optim$bestval,
  " m\n\n",
  sep = ""
)

calibration_predictions <- predict_transfer(
  large_data,
  calibration_rf,
  calibration_Kscale
)

large_data_with_prediction <- large_data %>%
  mutate(Z_MAX_transfer = calibration_predictions)

calibration_metrics <- calculate_metrics_cv(
  large_data_with_prediction,
  "Z_MAX_transfer"
)

validation_metrics <- calculate_metrics_cv(
  small_data_with_uncertainty,
  "Z_MAX_mean"
)

calibration_by_region <- large_data_with_prediction %>%
  group_by(Region) %>%
  group_modify(
    ~ calculate_metrics_cv(.x, "Z_MAX_transfer")
  ) %>%
  ungroup()

validation_by_region <- small_data_with_uncertainty %>%
  group_by(Region) %>%
  group_modify(
    ~ calculate_metrics_cv(.x, "Z_MAX_mean")
  ) %>%
  ungroup()

# ------------------------------------------------------------------------------
# 19. Print performance results
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("CALIBRATION PERFORMANCE\n")
cat("============================================================\n\n")
print(calibration_metrics)

cat("\n============================================================\n")
cat("VALIDATION PERFORMANCE\n")
cat("============================================================\n\n")
print(validation_metrics)

cat("\nCalibration by region:\n")
print(calibration_by_region)

cat("\nValidation by region:\n")
print(validation_by_region)

# ------------------------------------------------------------------------------
# 20. Save validation results
# ------------------------------------------------------------------------------

write.csv(
  small_data_with_uncertainty %>%
    select(
      Island,
      Region,
      Z_MAX_O,
      Z_MAX_mean,
      Z_MAX_median,
      Z_MAX_sd,
      Z_MAX_lower_pi_95,
      Z_MAX_upper_pi_95,
      Z_MAX_pi_width_95,
      covered_95,
      Z_MAX_error
    ),
  "transfer_validation_predictions_with_uncertainty.csv",
  row.names = FALSE
)

write.csv(
  coverage_overall,
  "transfer_validation_95PI_coverage_overall.csv",
  row.names = FALSE
)

write.csv(
  coverage_by_region,
  "transfer_validation_95PI_coverage_by_region.csv",
  row.names = FALSE
)

write.csv(
  region_counts,
  "transfer_validation_region_counts.csv",
  row.names = FALSE
)

performance_comparison <- bind_rows(
  calibration_metrics %>%
    mutate(Dataset = "Calibration (4 largest)"),
  validation_metrics %>%
    mutate(Dataset = "Validation (4 smallest)")
) %>%
  select(Dataset, everything())

write.csv(
  performance_comparison,
  "transfer_validation_performance_bootstrap.csv",
  row.names = FALSE
)

write.csv(
  calibration_by_region %>%
    mutate(Dataset = "Calibration"),
  "transfer_validation_calibration_by_region.csv",
  row.names = FALSE
)

write.csv(
  validation_by_region %>%
    mutate(Dataset = "Validation"),
  "transfer_validation_validation_by_region.csv",
  row.names = FALSE
)

cat("Saved validation prediction results.\n")
cat("Saved overall 95% PI coverage.\n")
cat("Saved regional 95% PI coverage.\n\n")

# ------------------------------------------------------------------------------
# 21. Bootstrap parameter summary
# ------------------------------------------------------------------------------


cat("\n============================================================\n")
cat("BOOTSTRAP PARAMETER SUMMARY\n")
cat("All 10,000 successful replicates\n")
cat("============================================================\n\n")

param_summary <- data.frame(
  Parameter = c(
    "RF (global)",
    "K scalar (global)"
  ),
  Mean = c(
    mean(boot_params_final$rf_global),
    mean(boot_params_final$Kscale_global)
  ),
  Median = c(
    median(boot_params_final$rf_global),
    median(boot_params_final$Kscale_global)
  ),
  SD = c(
    sd(boot_params_final$rf_global),
    sd(boot_params_final$Kscale_global)
  ),
  Lower_95 = c(
    quantile(
      boot_params_final$rf_global,
      0.025
    ),
    quantile(
      boot_params_final$Kscale_global,
      0.025
    )
  ),
  Upper_95 = c(
    quantile(
      boot_params_final$rf_global,
      0.975
    ),
    quantile(
      boot_params_final$Kscale_global,
      0.975
    )
  )
)

print(param_summary)

write.csv(
  param_summary,
  "bootstrap_Transfer_parameter_summary.csv",
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 22. Summary
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("TRANSFER VALIDATION WITH BOOTSTRAP COMPLETE\n")
cat("============================================================\n\n")

cat(
  "Calibration RMSE:   ",
  round(calibration_metrics$RMSE, 4),
  " m\n",
  sep = ""
)

cat(
  "Validation RMSE:    ",
  round(validation_metrics$RMSE, 4),
  " m\n",
  sep = ""
)

cat(
  "Validation Bias:    ",
  round(validation_metrics$Bias, 4),
  " m\n",
  sep = ""
)

cat(
  "Validation R²:      ",
  round(validation_metrics$R2, 4),
  "\n",
  sep = ""
)

cat(
  "\n95% PI coverage:    ",
  coverage_overall$n_covered,
  " / ",
  coverage_overall$n,
  " (",
  round(100 * coverage_overall$coverage, 1),
  "%)\n",
  sep = ""
)

cat(
  "Mean 95% PI width:  ",
  round(coverage_overall$mean_width, 3),
  " m\n",
  sep = ""
)

cat(
  "Median 95% PI width:",
  round(coverage_overall$median_width, 3),
  " m\n\n",
  sep = ""
)

cat("Files saved to Data_Files/:\n")
cat("  - bootstrap_Transfer_parameters.csv\n")
cat("  - bootstrap_Transfer_parameters_with_replicate_id.csv\n")
cat("  - bootstrap_Transfer_parameter_summary.csv\n")
cat("  - transfer_objective_verification.csv\n")
cat("  - transfer_validation_predictions_with_uncertainty.csv\n")
cat("  - transfer_validation_95PI_coverage_overall.csv\n")
cat("  - transfer_validation_95PI_coverage_by_region.csv\n")
cat("  - transfer_validation_region_counts.csv\n")
cat("  - transfer_validation_performance_bootstrap.csv\n")
cat("  - transfer_validation_calibration_by_region.csv\n")
cat("  - transfer_validation_validation_by_region.csv\n\n")

# ------------------------------------------------------------------------------
# 23. Plots
# ------------------------------------------------------------------------------

cat("Creating plots...\n\n")

plot_data <- small_data_with_uncertainty %>%
  arrange(Region, Island)

# Trace plots
cat("Creating trace plots for global parameters...\n")

windows()
par(mfrow = c(2, 1))
par(mar = c(3, 3, 2, 1), mgp = c(2, 0.7, 0))

plot(
  boot_params_with_id$rf_global,
  type = "l",
  xlab = "Bootstrap replicate",
  ylab = "Global Recharge Fraction",
  main = "Trace Plot: Global RF (10,000 bootstrap replicates)",
  las = 1,
  lwd = 0.5,
  col = "steelblue"
)

plot(
  boot_params_with_id$Kscale_global,
  type = "l",
  xlab = "Bootstrap replicate",
  ylab = "Global K Scalar",
  main = "Trace Plot: Global K Scalar (10,000 bootstrap replicates)",
  las = 1,
  lwd = 0.5,
  col = "darkred"
)

par(mfrow = c(1, 1))

# Parameter histograms
cat("Creating histograms for global parameters...\n")

windows()
par(mfrow = c(1, 2))
par(mar = c(3, 3, 2, 1), mgp = c(2, 0.7, 0))

hist(
  boot_params_mat$rf_global,
  freq = FALSE,
  xlim = c(0, 1),
  xlab = "Global Recharge Fraction",
  ylab = "Probability Density",
  main = "Distribution: Global RF\n(all 10,000 replicates)",
  las = 1,
  cex.main = 1,
  breaks = 40,
  col = "steelblue",
  border = "gray"
)

abline(
  v = mean(boot_params_mat$rf_global),
  col = "red",
  lwd = 2,
  lty = 2
)

abline(
  v = quantile(
    boot_params_mat$rf_global,
    c(0.025, 0.975)
  ),
  col = "orange",
  lwd = 1.5,
  lty = 3
)

legend(
  "topright",
  legend = c("Mean", "95% CI bounds"),
  col = c("red", "orange"),
  lty = c(2, 3),
  cex = 0.9
)

hist(
  boot_params_mat$Kscale_global,
  freq = FALSE,
  xlim = c(0, 4),
  xlab = "Global K Scalar",
  ylab = "Probability Density",
  main = "Distribution: Global K Scalar\n(all 10,000 replicates)",
  las = 1,
  cex.main = 1,
  breaks = 40,
  col = "darkred",
  border = "gray"
)

abline(
  v = mean(boot_params_mat$Kscale_global),
  col = "blue",
  lwd = 2,
  lty = 2
)

abline(
  v = quantile(
    boot_params_mat$Kscale_global,
    c(0.025, 0.975)
  ),
  col = "orange",
  lwd = 1.5,
  lty = 3
)

legend(
  "topright",
  legend = c("Mean", "95% CI bounds"),
  col = c("blue", "orange"),
  lty = c(2, 3),
  cex = 0.9
)

par(mfrow = c(1, 1))

# ==============================================================================
# Observed vs predicted validation plot
# Validation islands only, with 95% bootstrap prediction intervals
# ==============================================================================

cat("Creating observed vs predicted validation plot...\n")

windows()

# Validation region colors
region_levels <- unique(plot_data$Region)

color_map <- setNames(
  rainbow(length(region_levels)),
  region_levels
)

# Determine plotting range from the validation data
plot_max <- max(
  c(
    plot_data$Z_MAX_O,
    plot_data$Z_MAX_mean,
    plot_data$Z_MAX_lower_pi_95,
    plot_data$Z_MAX_upper_pi_95
  ),
  na.rm = TRUE
)

# Round up to a clean axis limit
plot_max <- ceiling(plot_max / 5) * 5

plot(
  NA,
  xlim = c(0, plot_max),
  ylim = c(0, plot_max),
  xlab = "Observed freshwater-lens depth (m)",
  ylab = "Predicted freshwater-lens depth (m)",
  main = paste(
    "Transfer Validation: Observed vs Predicted\n",
    "(95% bootstrap prediction intervals, n = 10,000)"
  ),
  las = 1
)

# 1:1 line
abline(
  a = 0,
  b = 1,
  lty = 2,
  col = "black",
  lwd = 1.5
)

# Add uncertainty intervals and validation points
for (reg in region_levels) {
  
  data_reg <- plot_data %>%
    filter(Region == reg)
  
  # 95% bootstrap prediction intervals
  for (i in seq_len(nrow(data_reg))) {
    
    segments(
      x0 = data_reg$Z_MAX_O[i],
      y0 = data_reg$Z_MAX_lower_pi_95[i],
      x1 = data_reg$Z_MAX_O[i],
      y1 = data_reg$Z_MAX_upper_pi_95[i],
      col = adjustcolor(
        color_map[reg],
        alpha.f = 0.5
      ),
      lwd = 1.5
    )
  }
  
  # Validation islands
  points(
    data_reg$Z_MAX_O,
    data_reg$Z_MAX_mean,
    pch = 17,
    col = color_map[reg],
    cex = 1.5
  )
}

# Region legend
legend(
  "topleft",
  legend = region_levels,
  col = color_map[region_levels],
  pch = 17,
  cex = 1,
  bty = "n"
)

# Prediction error distribution
cat("Creating prediction error distribution plot...\n")

windows()

hist(
  plot_data$Z_MAX_error,
  main = "Distribution of Prediction Errors\n(Mean - Observed)",
  xlab = "Prediction error (m)",
  ylab = "Frequency",
  breaks = 15,
  col = "lightblue",
  border = "gray"
)

abline(
  v = 0,
  col = "red",
  lwd = 2,
  lty = 2
)

# Uncertainty by validation island
cat("Creating uncertainty by island plot...\n")

windows()

barplot(
  plot_data$Z_MAX_sd,
  names.arg = paste(
    plot_data$Island,
    plot_data$Region,
    sep = "-"
  ),
  main = "Bootstrap Prediction Uncertainty (1 SD)\nby Validation Island",
  ylab = "Standard deviation (m)",
  las = 2,
  cex.names = 0.7
)

cat("All plots created.\n\n")

# ------------------------------------------------------------------------------
# 24. Final message
# ------------------------------------------------------------------------------

cat("============================================================\n")
cat("ALL ANALYSIS COMPLETE!\n")
cat("============================================================\n\n")

cat(
  "Elapsed bootstrap time: ",
  round(as.numeric(elapsed_time), 2),
  " hours\n",
  sep = ""
)

cat("Objective-equivalence test: PASSED\n")
cat("Bootstrap replicates: ", n_boot, "\n", sep = "")
cat(
  "Successful bootstrap replicates: ",
  n_successful,
  "\n",
  sep = ""
)

cat(
  "Validation prediction replicates: all ",
  n_final_boot,
  " successful bootstrap replicates\n",
  sep = ""
)

cat(
  "95% PI coverage: ",
  coverage_overall$n_covered,
  " / ",
  coverage_overall$n,
  " (",
  round(100 * coverage_overall$coverage, 1),
  "%)\n\n",
  sep = ""
)

cat("Plots generated:\n")
cat("  1. Trace plots: RF and K scalar over bootstrap replicates\n")
cat("  2. Histograms: RF and K scalar\n")
cat("  3. Observed vs predicted: validation islands with 95% prediction intervals\n")
cat("  4. Prediction errors\n")
cat("  5. Prediction uncertainty by validation island\n")

cat("\n============================================================\n")
cat("END OF TRANSFER VALIDATION\n")
cat("============================================================\n")



























# ==============================================================================
# Publication-quality transfer validation plot using ggplot2
# ==============================================================================
cat("Creating ggplot transfer validation figure...\n")

library(ggplot2)

# Set region order
plot_data$Region <- factor(
  plot_data$Region,
  levels = c(
    "Central_Indian",
    "Cook_Islands",
    "East_Indian",
    "North_Pacific"
  )
)

# Determine common axis limit
plot_max <- max(
  c(
    plot_data$Z_MAX_O,
    plot_data$Z_MAX_mean,
    plot_data$Z_MAX_lower_pi_95,
    plot_data$Z_MAX_upper_pi_95
  ),
  na.rm = TRUE
)

# Round up to the next 5 m and add 5 m of headroom
plot_max <- ceiling(plot_max / 5) * 5 + 5

# Create plot
transfer_plot <- ggplot(
  plot_data,
  aes(
    x = Z_MAX_O,
    y = Z_MAX_mean,
    color = Region
  )
) +
  
  # 1:1 line
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.6,
    color = "gray40"
  ) +
  
  # 95% bootstrap prediction intervals
  geom_errorbar(
    aes(
      ymin = Z_MAX_lower_pi_95,
      ymax = Z_MAX_upper_pi_95
    ),
    width = 0,
    linewidth = 0.6,
    alpha = 0.65
  ) +
  
  # Validation observations
  geom_point(
    size = 3,
    shape = 16
  ) +
  
  # Equal scales with extended axes
  coord_equal(
    xlim = c(0, plot_max),
    ylim = c(0, plot_max),
    expand = FALSE
  ) +
  
  # Region colors
  scale_color_brewer(
    palette = "Dark2"
  ) +
  
  # Labels
  labs(
    x = "Observed freshwater-lens depth (m)",
    y = "Predicted freshwater-lens depth (m)",
    color = "Region"
  ) +
  
  # Clean theme
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5)
  )

print(transfer_plot)

# Save high-resolution version
ggsave(
  filename = "transfer_validation_ggplot.png",
  plot = transfer_plot,
  width = 15,
  height = 12,
  units = "cm",
  dpi = 600,
  bg = "white"
)

cat("ggplot transfer validation figure saved.\n")

