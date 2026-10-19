# ==============================================================================
# PRECIPITATION PRODUCT vs OBSERVED PRECIPITATION
# ==============================================================================

library(dplyr)
library(tidyr)
library(lubridate)

setwd("C:/Users/Owner/OneDrive/Desktop/Github")


# ==============================================================================
# 1. GRIDDED PRECIPITATION DATA
# ==============================================================================

CDR_df <- as.data.frame(CDR_monthly)
MSWEP_df <- as.data.frame(MSWEP_results_df)
NASA_df <- as.data.frame(NASA_monthly)
JAXA_df <- as.data.frame(JAXA_results_df)
ERA5_df <- as.data.frame(ERA5_monthly)

names(CDR_df) <- tolower(names(CDR_df))
names(MSWEP_df) <- tolower(names(MSWEP_df))
names(NASA_df) <- tolower(names(NASA_df))
names(JAXA_df) <- tolower(names(JAXA_df))
names(ERA5_df) <- tolower(names(ERA5_df))

CDR_df$island <- tolower(CDR_df$island)
MSWEP_df$island <- tolower(MSWEP_df$island)
NASA_df$island <- tolower(NASA_df$island)
JAXA_df$island <- tolower(JAXA_df$island)
ERA5_df$island <- tolower(ERA5_df$island)


# ==============================================================================
# 2. STANDARDIZE PRECIPITATION COLUMN
# ==============================================================================

find_precip_column <- function(df) {
  
  possible_names <- c(
    "monthly_precip",
    "monthly_precipitation",
    "precipitation",
    "precip",
    "precip_mm"
  )
  
  found <- possible_names[possible_names %in% names(df)]
  
  if (length(found) == 0) {
    stop(
      "Could not identify monthly precipitation column. Available columns: ",
      paste(names(df), collapse = ", ")
    )
  }
  
  found[1]
}


names(CDR_df)[
  names(CDR_df) == find_precip_column(CDR_df)
] <- "monthly_precip"

names(MSWEP_df)[
  names(MSWEP_df) == find_precip_column(MSWEP_df)
] <- "monthly_precip"

names(NASA_df)[
  names(NASA_df) == find_precip_column(NASA_df)
] <- "monthly_precip"

names(JAXA_df)[
  names(JAXA_df) == find_precip_column(JAXA_df)
] <- "monthly_precip"

names(ERA5_df)[
  names(ERA5_df) == find_precip_column(ERA5_df)
] <- "monthly_precip"


# ==============================================================================
# 3. LABEL PRODUCTS
# ==============================================================================

CDR_df$product <- "CDR"
MSWEP_df$product <- "MSWEP"
NASA_df$product <- "NASA"
JAXA_df$product <- "JAXA"
ERA5_df$product <- "ERA5"


gridded_data <- dplyr::bind_rows(
  CDR_df,
  MSWEP_df,
  NASA_df,
  JAXA_df,
  ERA5_df
)

gridded_data$year <- as.integer(gridded_data$year)


# ==============================================================================
# 4. CHECK GRIDDED DATA
# ==============================================================================

cat("\n")
cat("GRIDDED DATA SUMMARY\n")
cat("====================\n")

print(
  gridded_data %>%
    dplyr::group_by(product) %>%
    dplyr::summarise(
      n_rows = dplyr::n(),
      n_islands = dplyr::n_distinct(island),
      n_missing = sum(is.na(monthly_precip)),
      first_year = min(year, na.rm = TRUE),
      last_year = max(year, na.rm = TRUE),
      mean_monthly_precip = mean(
        monthly_precip,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
)


# ==============================================================================
# 5. READ LOCAL OBSERVATIONS
# ==============================================================================

laura <- read.csv("Data_Files/Local_Rainfall/Laura.csv")

cocos <- read.csv("Data_Files/Local_Rainfall/Cocos.csv")

diego <- read.csv("Data_Files/Local_Rainfall/Diego_Garcia.csv")

ulithi <- read.csv("Data_Files/Local_Rainfall/Ulithi.csv")

tarawa <- read.csv("Data_Files/Local_Rainfall/Tarawa.csv")

#MALDIVES---------------

gan_addu <- read.csv("Data_Files/Local_Rainfall/Gan_Addu.csv", skip = 4)

hanimaadhoo <- read.csv("Data_Files/Local_Rainfall/Hanimaadhoo.csv", skip = 4)

kadhdhoo <- read.csv("Data_Files/Local_Rainfall/Kadhdhoo.csv", skip = 4)


# ==============================================================================
# 6. RAW DATA CHECK
# ==============================================================================

cat("\n")
cat("RAW OBSERVED PRECIPITATION\n")
cat("==========================\n")

cat("\nLaura:\n")
print(summary(laura$PRCP))

cat("\nCocos:\n")
print(summary(cocos$PRCP))

cat("\nDiego Garcia:\n")
print(summary(diego$PRCP))

cat("\nUlithi/Pingelap:\n")
print(summary(ulithi$PRCP))

cat("\nTarawa:\n")
print(summary(tarawa$PRCP))


# ==============================================================================
# 7. DAILY -> MONTHLY OBSERVED PRECIPITATION
# ==============================================================================
#
# IMPORTANT:
#
# PRCP is ALREADY IN INCHES.
#
# Therefore:
#
#   monthly_inches = sum(daily PRCP)
#
#   monthly_mm = monthly_inches * 25.4
#
# A month is considered valid when >=90% of its daily records have
# non-missing precipitation.
#
# ==============================================================================


daily_to_monthly <- function(
    df,
    island_name,
    station_id = NULL
) {
  
  # Filter to requested station
  if (
    !is.null(station_id) &&
    "STATION" %in% names(df)
  ) {
    df <- df[
      df$STATION == station_id,
      ,
      drop = FALSE
    ]
  }
  
  # Convert dates
  df$DATE <- as.Date(as.character(df$DATE))
  
  df <- df[
    !is.na(df$DATE),
    ,
    drop = FALSE
  ]
  
  # Extract year and month
  df$year <- as.integer(format(df$DATE, "%Y"))
  df$month <- as.integer(format(df$DATE, "%m"))
  
  # PRCP is reported in inches
  df$PRCP_inches <- as.numeric(df$PRCP)
  
  # Aggregate to monthly totals.
  # Require at least 90% of calendar days to have valid PRCP.
  monthly <- df %>%
    dplyr::group_by(year, month) %>%
    dplyr::summarise(
      n_daily_records = dplyr::n(),
      
      n_valid_rainfall = sum(
        !is.na(PRCP_inches)
      ),
      
      days_expected = lubridate::days_in_month(
        lubridate::make_date(
          year[1],
          month[1],
          1
        )
      ),
      
      monthly_precip_inches = ifelse(
        n_valid_rainfall / days_expected >= 0.90,
        sum(PRCP_inches, na.rm = TRUE),
        NA_real_
      ),
      
      .groups = "drop"
    )
  
  # Convert inches to millimeters
  monthly$monthly_precip <-
    monthly$monthly_precip_inches * 25.4
  
  # Add island identifier
  monthly$island <- tolower(island_name)
  
  return(monthly)
}


# ==============================================================================
# 8. BUILD OBSERVED DATA
# ==============================================================================

observed_data <- dplyr::bind_rows(
  
  daily_to_monthly(
    laura,
    "laura",
    "RMC00914399"
  ),
  
  daily_to_monthly(
    cocos,
    "cocos",
    "CK000096996"
  ),
  
  daily_to_monthly(
    diego,
    "cantonment",
    "IOW00070701"
  ),
  
  daily_to_monthly(
    ulithi,
    "falalop",
    "FMC00914892"
  ),
  
  daily_to_monthly(
    tarawa,
    "tarawa",
    "KR000091610"
  )
)


# ==============================================================================
# 8B. MALDIVES DAILY -> MONTHLY OBSERVED PRECIPITATION
# ==============================================================================
#
# Maldives station files report:
#
#   Daily Accumulated Rainfall (mm)
#
# Therefore:
#
#   monthly_precip_mm = sum(daily rainfall)
#
# No inches -> mm conversion is needed.
#
# A month is considered valid when >=90% of its daily records
# contain non-missing rainfall.
#
# ==============================================================================


maldives_daily_to_monthly <- function(
    df,
    island_name
) {
  
  # Standardize column names
  names(df) <- tolower(trimws(names(df)))
  
  # Identify date and rainfall columns
  date_col <- names(df)[
    grepl("date", names(df))
  ][1]
  
  precip_col <- names(df)[
    grepl("rain|precip", names(df))
  ][1]
  
  if (
    length(date_col) == 0 ||
    is.na(date_col) ||
    length(precip_col) == 0 ||
    is.na(precip_col)
  ) {
    stop(
      "Could not identify date/rainfall columns for ",
      island_name,
      ". Available columns: ",
      paste(names(df), collapse = ", ")
    )
  }
  
  # Preserve original values for quality checks
  date_character <- trimws(
    as.character(df[[date_col]])
  )
  
  rainfall_character <- trimws(
    as.character(df[[precip_col]])
  )
  
  # Standardize slash separators.
  # Dates are interpreted as day-month-year.
  date_character <- gsub(
    "/",
    "-",
    date_character
  )
  
  parsed_date <- as.Date(
    date_character,
    format = "%d-%m-%Y"
  )
  
  # Report unparseable dates before dropping them
  n_bad_dates <- sum(is.na(parsed_date))
  
  cat(
    "\nMaldives station:", island_name,
    "\nRows in source file:", nrow(df),
    "\nUnparseable dates:", n_bad_dates,
    "\n"
  )
  
  # Attach parsed dates and rainfall values
  df$DATE <- parsed_date
  
  # Normalize trace labels
  trace_values <- toupper(
    rainfall_character
  ) %in% c("TR", "TRACE")
  
  # Convert rainfall to numeric.
  # Explicitly treat trace rainfall as zero for the
  # primary analysis, rather than as missing.
  rainfall_character[
    trace_values
  ] <- "0"
  
  df$rainfall_mm <- suppressWarnings(
    as.numeric(rainfall_character)
  )
  
  # Treat blank strings as missing
  df$rainfall_mm[
    rainfall_character == ""
  ] <- NA_real_
  
  # Report trace and unparseable rainfall values
  n_trace <- sum(trace_values, na.rm = TRUE)
  
  n_bad_rainfall <- sum(
    is.na(df$rainfall_mm) &
      !is.na(rainfall_character) &
      rainfall_character != ""
  )
  
  cat(
    "Trace rainfall records treated as 0 mm:",
    n_trace,
    "\nUnrecognized rainfall values:",
    n_bad_rainfall,
    "\n"
  )
  
  # Remove rows with invalid dates
  df <- df %>%
    dplyr::filter(!is.na(DATE))
  
  # Extract year and month
  df$year <- as.integer(
    format(df$DATE, "%Y")
  )
  
  df$month <- as.integer(
    format(df$DATE, "%m")
  )
  
  # Check for duplicate dates before aggregation
  duplicate_dates <- df %>%
    dplyr::count(DATE) %>%
    dplyr::filter(n > 1)
  
  if (nrow(duplicate_dates) > 0) {
    warning(
      "Duplicate dates found for ",
      island_name,
      ". Check the source file before interpreting totals."
    )
  }
  
  # Aggregate to monthly totals.
  #
  # A month is considered valid when at least 90% of the
  # calendar days have valid rainfall observations.
  #
  # "TR" values have already been converted to 0 mm above
  # and therefore count as valid observations.
  
  monthly <- df %>%
    dplyr::group_by(year, month) %>%
    dplyr::summarise(
      
      # Number of records in the source file
      n_daily_records = dplyr::n(),
      
      # Number of distinct calendar dates represented
      n_dates_present = dplyr::n_distinct(DATE),
      
      # Number of distinct calendar dates with valid rainfall
      n_valid_rainfall = dplyr::n_distinct(
        DATE[!is.na(rainfall_mm)]
      ),
      
      # Number of calendar days expected in the month
      days_expected = lubridate::days_in_month(
        lubridate::make_date(
          year[1],
          month[1],
          1
        )
      ),
      
      # Monthly total is retained only if >=90% of
      # calendar days have valid observations
      monthly_precip = ifelse(
        n_valid_rainfall / days_expected >= 0.90,
        sum(rainfall_mm, na.rm = TRUE),
        NA_real_
      ),
      
      .groups = "drop"
    )
  
  # Add station/island identifier
  monthly$island <- tolower(island_name)
  
  return(monthly)
}


# ==============================================================================
# 8C. BUILD MALDIVES OBSERVED DATA
# ==============================================================================

maldives_observed <- dplyr::bind_rows(
  
  maldives_daily_to_monthly(
    gan_addu,
    "gan (addu)"
  ),
  
  maldives_daily_to_monthly(
    hanimaadhoo,
    "hanimaadhoo"
  ),
  
  maldives_daily_to_monthly(
    kadhdhoo,
    "gan (laamu)"
  )
)


# Add Maldives station observations
observed_data <- dplyr::bind_rows(
  observed_data,
  maldives_observed
)

# ==============================================================================
# 9. IMMEDIATELY CHECK OBSERVED MONTHLY DATA
# ==============================================================================

cat("\n")
cat("OBSERVED MONTHLY DATA CHECK\n")
cat("===========================\n")

print(
  observed_data %>%
    dplyr::group_by(island) %>%
    dplyr::summarise(
      
      n_months = dplyr::n(),
      
      n_valid_months =
        sum(!is.na(monthly_precip)),
      
      first_year =
        min(year),
      
      last_year =
        max(year),
      
      mean_monthly_inches =
        mean(
          monthly_precip_inches,
          na.rm = TRUE
        ),
      
      mean_monthly_mm =
        mean(
          monthly_precip,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    )
)

# ==============================================================================
# 9B. MALDIVES OBSERVED MONTHLY DATA CHECK
# ==============================================================================

cat("\n")
cat("MALDIVES OBSERVED MONTHLY DATA CHECK\n")
cat("====================================\n")

print(
  maldives_observed %>%
    dplyr::group_by(island) %>%
    dplyr::summarise(
      
      n_months =
        dplyr::n(),
      
      n_valid_months =
        sum(!is.na(monthly_precip)),
      
      first_year =
        min(year),
      
      last_year =
        max(year),
      
      mean_monthly_mm =
        mean(
          monthly_precip,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    )
)


# ==============================================================================
# 9C. MALDIVES MONTHLY COMPLETENESS BY YEAR
# ==============================================================================

cat("\n")
cat("MALDIVES MONTHLY COMPLETENESS BY YEAR\n")
cat("=====================================\n")

maldives_year_validity <- maldives_observed %>%
  dplyr::group_by(island, year) %>%
  dplyr::summarise(
    
    n_months_present = dplyr::n_distinct(month),
    
    n_valid_months = sum(
      !is.na(monthly_precip)
    ),
    
    missing_or_invalid_months = 12 - n_valid_months,
    
    valid_year = (
      n_months_present == 12 &
        n_valid_months == 12
    ),
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(island, year)

print(maldives_year_validity, n = Inf)

cat("\nNUMBER OF VALID YEARS BY MALDIVES STATION\n")
cat("=========================================\n")

print(
  maldives_year_validity %>%
    dplyr::group_by(island) %>%
    dplyr::summarise(
      n_years_in_source = dplyr::n(),
      
      n_valid_years = sum(valid_year),
      
      first_valid_year = ifelse(
        any(valid_year),
        min(year[valid_year]),
        NA_integer_
      ),
      
      last_valid_year = ifelse(
        any(valid_year),
        max(year[valid_year]),
        NA_integer_
      ),
      
      .groups = "drop"
    )
)

# Save the complete year-by-year diagnostic
write.csv(
  maldives_year_validity,
  "Data_Files/maldives_year_validity.csv",
  row.names = FALSE
)

# ==============================================================================
# 10. OBSERVED MONTHLY -> ANNUAL
# ==============================================================================
#
# Require all 12 months to be valid.
#
# ==============================================================================


observed_annual <- observed_data %>%
  
  dplyr::group_by(island, year) %>%
  
  dplyr::summarise(
    
    n_months_present = dplyr::n_distinct(month),
    
    n_valid_months = sum(
      !is.na(monthly_precip)
    ),
    
    annual_precip_obs = ifelse(
      n_months_present == 12 &
        n_valid_months == 12,
      
      sum(monthly_precip, na.rm = TRUE),
      
      NA_real_
    ),
    
    .groups = "drop"
  ) %>%
  
  dplyr::filter(
    !is.na(annual_precip_obs)
  )

# ==============================================================================
# 11. CRITICAL OBSERVED ANNUAL CHECK
# ==============================================================================

cat("\n")
cat("OBSERVED ANNUAL PRECIPITATION CHECK\n")
cat("===================================\n")

observed_summary <- observed_annual %>%
  
  dplyr::group_by(island) %>%
  
  dplyr::summarise(
    
    n_years =
      dplyr::n(),
    
    first_year =
      min(year),
    
    last_year =
      max(year),
    
    mean_annual_mm =
      mean(annual_precip_obs),
    
    min_annual_mm =
      min(annual_precip_obs),
    
    max_annual_mm =
      max(annual_precip_obs),
    
    .groups = "drop"
  )

print(observed_summary)


# ==============================================================================
# 12. TARGET ISLANDS / YEARS
# ==============================================================================

target_islands <- c(
  "laura",
  "cantonment",
  "falalop",
  "gan (addu)",
  "hanimaadhoo",
  "gan (laamu)"
)

target_years <- 1998:2020


gridded_filt <- gridded_data %>%
  
  dplyr::filter(
    island %in% target_islands,
    year %in% target_years
  )


observed_filt <- observed_annual %>%
  
  dplyr::filter(
    island %in% target_islands,
    year %in% target_years
  )


# ==============================================================================
# 12B. VALID OBSERVED YEARS USED IN THE COMPARISON
# ==============================================================================

cat("\n")
cat("VALID OBSERVED YEARS USED IN THE COMPARISON\n")
cat("============================================\n")

print(
  observed_filt %>%
    dplyr::filter(
      island %in% c(
        "gan (addu)",
        "hanimaadhoo",
        "gan (laamu)"
      )
    ) %>%
    dplyr::group_by(island) %>%
    dplyr::summarise(
      n_years = dplyr::n(),
      years_included = paste(
        sort(year),
        collapse = ", "
      ),
      .groups = "drop"
    )
)

# ==============================================================================
# 13. GRIDDED MONTHLY -> ANNUAL
# ==============================================================================

gridded_annual <- gridded_filt %>%
  
  dplyr::group_by(
    island,
    product,
    year
  ) %>%
  
  dplyr::summarise(
    
    n_months =
      dplyr::n_distinct(month),
    
    n_valid_months =
      sum(!is.na(monthly_precip)),
    
    annual_precip =
      ifelse(
        n_months == 12 &
          n_valid_months == 12,
        
        sum(
          monthly_precip,
          na.rm = TRUE
        ),
        
        NA_real_
      ),
    
    .groups = "drop"
  ) %>%
  
  dplyr::filter(
    !is.na(annual_precip)
  )


# ==============================================================================
# 14. JOIN OBSERVATIONS AND GRIDDED PRODUCTS
# ==============================================================================

individual_comparison <- gridded_annual %>%
  
  dplyr::inner_join(
    observed_filt,
    by = c(
      "island",
      "year"
    )
  )


# ==============================================================================
# 15. INDIVIDUAL PRODUCT PERFORMANCE
# ==============================================================================

products <- c(
  "CDR",
  "MSWEP",
  "NASA",
  "JAXA",
  "ERA5"
)


product_comparisons <- list()

product_stats <- data.frame()


for (prod in products) {
  
  temp <- individual_comparison %>%
    
    dplyr::filter(
      product == prod,
      !is.na(annual_precip),
      !is.na(annual_precip_obs)
    )
  
  
  temp$precip_error <-
    temp$annual_precip -
    temp$annual_precip_obs
  
  
  product_comparisons[[prod]] <- temp
  
  
  if (nrow(temp) > 0) {
    
    obs_mean <-
      mean(temp$annual_precip_obs)
    
    product_mean <-
      mean(temp$annual_precip)
    
    bias <-
      mean(temp$precip_error)
    
    MAE <-
      mean(abs(temp$precip_error))
    
    RMSE <-
      sqrt(mean(temp$precip_error^2))
    
    correlation <-
      cor(
        temp$annual_precip,
        temp$annual_precip_obs,
        use = "complete.obs"
      )
    
    
    product_stats <- dplyr::bind_rows(
      
      product_stats,
      
      data.frame(
        
        product = prod,
        
        n_island_years = nrow(temp),
        
        obs_mean = obs_mean,
        
        product_mean = product_mean,
        
        bias = bias,
        
        percent_bias =
          (bias / obs_mean) * 100,
        
        MAE = MAE,
        
        RMSE = RMSE,
        
        correlation = correlation
      )
    )
  }
}


# ==============================================================================
# 16. FIVE-PRODUCT ENSEMBLE
# ==============================================================================

gridded_wide <- gridded_annual %>%
  
  dplyr::select(
    island,
    year,
    product,
    annual_precip
  ) %>%
  
  tidyr::pivot_wider(
    id_cols = c(
      island,
      year
    ),
    names_from = product,
    values_from = annual_precip
  )


ensemble_comparison <- gridded_wide %>%
  
  dplyr::inner_join(
    observed_filt,
    by = c(
      "island",
      "year"
    )
  ) %>%
  
  dplyr::filter(
    !is.na(CDR),
    !is.na(MSWEP),
    !is.na(NASA),
    !is.na(JAXA),
    !is.na(ERA5),
    !is.na(annual_precip_obs)
  )


ensemble_comparison$ensemble_mean <-
  rowMeans(
    ensemble_comparison[
      ,
      c(
        "CDR",
        "MSWEP",
        "NASA",
        "JAXA",
        "ERA5"
      )
    ]
  )


# ==============================================================================
# 17. OVERALL ENSEMBLE STATISTICS
# ==============================================================================

if (nrow(ensemble_comparison) > 0) {
  
  ensemble_errors <-
    ensemble_comparison$ensemble_mean -
    ensemble_comparison$annual_precip_obs
  
  obs_m <-
    mean(
      ensemble_comparison$annual_precip_obs
    )
  
  ens_m <-
    mean(
      ensemble_comparison$ensemble_mean
    )
  
  b <-
    mean(ensemble_errors)
  
  m <-
    mean(abs(ensemble_errors))
  
  r <-
    sqrt(mean(ensemble_errors^2))
  
  c <-
    cor(
      ensemble_comparison$ensemble_mean,
      ensemble_comparison$annual_precip_obs,
      use = "complete.obs"
    )
  
  
  overall_stats <- data.frame(
    
    n_island_years =
      nrow(ensemble_comparison),
    
    obs_mean =
      obs_m,
    
    ens_mean =
      ens_m,
    
    bias =
      b,
    
    percent_bias =
      (b / obs_m) * 100,
    
    MAE =
      m,
    
    RMSE =
      r,
    
    correlation =
      c
  )
  
} else {
  
  overall_stats <- data.frame(
    n_island_years = 0,
    obs_mean = NA_real_,
    ens_mean = NA_real_,
    bias = NA_real_,
    percent_bias = NA_real_,
    MAE = NA_real_,
    RMSE = NA_real_,
    correlation = NA_real_
  )
}


# ==============================================================================
# 18. ENSEMBLE STATISTICS BY ISLAND
# ==============================================================================

ensemble_by_island <- ensemble_comparison %>%
  
  dplyr::group_by(island) %>%
  
  dplyr::summarise(
    
    n_island_years =
      dplyr::n(),
    
    first_year =
      min(year),
    
    last_year =
      max(year),
    
    obs_mean =
      mean(annual_precip_obs),
    
    ensemble_mean_mm =
      mean(ensemble_mean),
    
    bias =
      mean(
        ensemble_mean -
          annual_precip_obs
      ),
    
    percent_bias =
      (
        mean(
          ensemble_mean -
            annual_precip_obs
        ) /
          mean(annual_precip_obs)
      ) * 100,
    
    MAE =
      mean(
        abs(
          ensemble_mean -
            annual_precip_obs
        )
      ),
    
    RMSE =
      sqrt(
        mean(
          (
            ensemble_mean -
              annual_precip_obs
          )^2
        )
      ),
    
    correlation =
      cor(
        ensemble_mean,
        annual_precip_obs,
        use = "complete.obs"
      ),
    
    .groups = "drop"
  )


# ==============================================================================
# 19. ENSEMBLE MEAN BY ISLAND
# ==============================================================================

ensemble_mean_comparison <- ensemble_comparison %>%
  
  dplyr::group_by(island) %>%
  
  dplyr::summarise(
    
    n_years =
      dplyr::n(),
    
    first_year =
      min(year),
    
    last_year =
      max(year),
    
    observed_mean_mm =
      mean(annual_precip_obs),
    
    ensemble_mean_mm =
      mean(ensemble_mean),
    
    bias_mm =
      ensemble_mean_mm -
      observed_mean_mm,
    
    percent_bias =
      (
        bias_mm /
          observed_mean_mm
      ) * 100,
    
    .groups = "drop"
  )


# ==============================================================================
# 20. COVERAGE
# ==============================================================================

product_coverage <- gridded_annual %>%
  
  dplyr::group_by(product) %>%
  
  dplyr::summarise(
    
    n_island_years =
      dplyr::n(),
    
    n_islands =
      dplyr::n_distinct(island),
    
    first_year =
      min(year),
    
    last_year =
      max(year),
    
    .groups = "drop"
  )


ensemble_coverage <- ensemble_comparison %>%
  
  dplyr::group_by(island) %>%
  
  dplyr::summarise(
    
    n_years =
      dplyr::n(),
    
    first_year =
      min(year),
    
    last_year =
      max(year),
    
    .groups = "drop"
  )


# ==============================================================================
# 21. SAVE OUTPUTS
# ==============================================================================

for (prod in products) {
  
  write.csv(
    product_comparisons[[prod]],
    paste0(
      "Data_Files/precipitation_comparison_",
      prod,
      ".csv"
    ),
    row.names = FALSE
  )
}


write.csv(
  individual_comparison,
  "Data_Files/precipitation_comparison_all_products.csv",
  row.names = FALSE
)

write.csv(
  ensemble_comparison,
  "Data_Files/precipitation_five_product_ensemble.csv",
  row.names = FALSE
)

write.csv(
  overall_stats,
  "Data_Files/overall_stats.csv",
  row.names = FALSE
)

write.csv(
  product_stats,
  "Data_Files/product_stats.csv",
  row.names = FALSE
)

write.csv(
  ensemble_by_island,
  "Data_Files/ensemble_stats_by_island.csv",
  row.names = FALSE
)

write.csv(
  ensemble_mean_comparison,
  "Data_Files/ensemble_mean_comparison.csv",
  row.names = FALSE
)

write.csv(
  product_coverage,
  "Data_Files/product_coverage.csv",
  row.names = FALSE
)

write.csv(
  ensemble_coverage,
  "Data_Files/ensemble_coverage.csv",
  row.names = FALSE
)


# ==============================================================================
# 22. FINAL RESULTS
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("OBSERVED ANNUAL PRECIPITATION\n")
cat("================================================================================\n")

print(observed_summary)


cat("\n")
cat("================================================================================\n")
cat("INDIVIDUAL PRODUCT PERFORMANCE\n")
cat("================================================================================\n")

print(product_stats)


cat("\n")
cat("================================================================================\n")
cat("FIVE-PRODUCT ENSEMBLE PERFORMANCE\n")
cat("================================================================================\n")

print(overall_stats)


cat("\n")
cat("================================================================================\n")
cat("ENSEMBLE BY ISLAND\n")
cat("================================================================================\n")

print(ensemble_by_island)


cat("\n")
cat("================================================================================\n")
cat("ENSEMBLE MEAN COMPARISON\n")
cat("================================================================================\n")

print(ensemble_mean_comparison)


cat("\n")
cat("================================================================================\n")
cat("PRODUCT COVERAGE\n")
cat("================================================================================\n")

print(product_coverage)


cat("\n")
cat("================================================================================\n")
cat("ENSEMBLE COVERAGE\n")
cat("================================================================================\n")

print(ensemble_coverage)


cat("\n")
cat("================================================================================\n")
cat("DONE\n")
cat("================================================================================\n")

