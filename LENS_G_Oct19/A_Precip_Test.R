# ============================================================
# PRECIPITATION EXTRACTION FOR ATOLL FRESHWATER-LENS MODEL
# ============================================================
#
# Products:
#   NASA
#   ERA5
#   JAXA
#   MSWEP
#   CDR
#
# Output:
#   Monthly precipitation by island and product
#   Annual precipitation by island and product
#
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

# Install once if necessary:
#
# install.packages("hdf5r")
# install.packages("ncdf4")
# install.packages("dplyr")
# install.packages("lubridate")
# install.packages("purrr")
# install.packages("tidyr")
# install.packages("readr")

library(hdf5r)
library(ncdf4)
library(dplyr)
library(lubridate)
library(purrr)
library(tidyr)
library(readr)


# ============================================================
# 2. WORKING DIRECTORY
# ============================================================

setwd(
  "C:/Users/Owner/OneDrive/Desktop/Github"
)


# ============================================================
# 3. DIRECTORIES
# ============================================================

precip_dir <-
  "C:/Users/Owner/OneDrive/Desktop/BAILEY/Precipitation"

rds_dir <-
  "Data_Files/RDS"

output_dir <-
  "Data_Files/Precipitation"

if (!dir.exists(rds_dir)) {
  
  dir.create(
    rds_dir,
    recursive = TRUE
  )
}

if (!dir.exists(output_dir)) {
  
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


# ============================================================
# 4. READ ISLAND COORDINATES
# ============================================================

final_results <- read.csv(
  "Data_Files/final_results.csv",
  stringsAsFactors = FALSE
)

island_features <-
  final_results %>%
  
  dplyr::filter(
    !is.na(Latitude),
    !is.na(Longitude),
    !is.na(Island)
  ) %>%
  
  dplyr::select(
    Island,
    Latitude,
    Longitude
  ) %>%
  
  dplyr::distinct()

cat(
  "Number of islands:",
  nrow(island_features),
  "\n"
)


# ============================================================
# 5. DATE RANGES
# ============================================================

start_date <-
  as.Date("1990-01-01")

end_date <-
  as.Date("2024-12-01")

satellite_start_date <-
  as.Date("1998-01-01")

mswep_end_date <-
  as.Date("2019-12-01")

cdr_end_date <-
  as.Date("2023-12-01")


# ============================================================
# 6. HELPER FUNCTIONS
# ============================================================


# ------------------------------------------------------------
# Find nearest latitude
# ------------------------------------------------------------

nearest_lat_index <- function(
    lat_values,
    target_lat
) {
  
  return(
    which.min(
      abs(
        lat_values -
          target_lat
      )
    )
  )
}


# ------------------------------------------------------------
# Adjust longitude depending on grid convention
#
# If grid is 0 to 360:
#   -162 becomes 198
#
# If grid is -180 to 180:
#   198 becomes -162
# ------------------------------------------------------------

adjust_longitude <- function(
    lon_values,
    target_lon
) {
  
  lon_min <-
    min(
      lon_values,
      na.rm = TRUE
    )
  
  lon_max <-
    max(
      lon_values,
      na.rm = TRUE
    )
  
  if (
    lon_min >= 0 &&
    lon_max > 180
  ) {
    
    target_lon <-
      target_lon %% 360
    
    if (
      target_lon < 0
    ) {
      
      target_lon <-
        target_lon + 360
    }
    
  } else {
    
    target_lon <-
      (
        (
          target_lon + 180
        ) %% 360
      ) - 180
  }
  
  return(
    target_lon
  )
}


# ------------------------------------------------------------
# Find nearest longitude
# ------------------------------------------------------------

nearest_lon_index <- function(
    lon_values,
    target_lon
) {
  
  target_lon_adjusted <-
    adjust_longitude(
      lon_values,
      target_lon
    )
  
  return(
    which.min(
      abs(
        lon_values -
          target_lon_adjusted
      )
    )
  )
}


# ------------------------------------------------------------
# Find NetCDF files
# ------------------------------------------------------------

find_nc_files <- function(
    directory
) {
  
  list.files(
    directory,
    pattern = "\\.nc$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
}


# ------------------------------------------------------------
# Find HDF5 files
# ------------------------------------------------------------

find_hdf5_files <- function(
    directory
) {
  
  list.files(
    directory,
    pattern = "\\.(h5|hdf5)$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
}


# ------------------------------------------------------------
# Extract first YYYYMMDD or YYYYMM date from filename
# ------------------------------------------------------------

extract_date_from_filename <- function(
    filename
) {
  
  basename_only <-
    basename(
      filename
    )
  
  # First look for YYYYMMDD
  match_8 <-
    regexpr(
      "[0-9]{8}",
      basename_only
    )
  
  if (
    match_8[1] != -1
  ) {
    
    date_string <-
      regmatches(
        basename_only,
        match_8
      )
    
    parsed_date <-
      as.Date(
        date_string,
        format = "%Y%m%d"
      )
    
    if (
      !is.na(parsed_date)
    ) {
      
      return(
        parsed_date
      )
    }
  }
  
  
  # Then look for YYYYMM
  match_6 <-
    regexpr(
      "[0-9]{6}",
      basename_only
    )
  
  if (
    match_6[1] != -1
  ) {
    
    date_string <-
      regmatches(
        basename_only,
        match_6
      )
    
    parsed_date <-
      as.Date(
        paste0(
          date_string,
          "01"
        ),
        format = "%Y%m%d"
      )
    
    if (
      !is.na(parsed_date)
    ) {
      
      return(
        parsed_date
      )
    }
  }
  
  
  return(
    as.Date(NA)
  )
}


# ============================================================
# 7. PRE-COMPUTE GRID INDICES
# ============================================================
#
# These are calculated separately for each product because
# each product has a different grid.
#


# ============================================================
# 8. NASA
# ============================================================


NASA_dir <- "C:/Users/Owner/OneDrive/Desktop/BAILEY/Precipitation/NASA"

NASA_files <- sort(list.files(
  NASA_dir,
  pattern = "\\.HDF5$",
  full.names = TRUE
))

# Dates from filenames
nasa_years <- as.Date(
  substr(basename(NASA_files), 1, 8),
  format = "%Y%m%d"
)

stopifnot(length(NASA_files) == length(nasa_years))


# ------------------------------------------------------------
# Find nearest NASA grid cell for each island ONCE
# ------------------------------------------------------------

# We can get the coordinate grids from the first file
h5f <- hdf5r::H5File$new(NASA_files[1], mode = "r")

NASA_lat <- h5f[["/Grid/lat"]]$read()
NASA_lon <- h5f[["/Grid/lon"]]$read()

h5f$close_all()


NASA_indices <- island_features %>%
  dplyr::mutate(
    lat_idx = purrr::map_int(
      Latitude,
      ~ which.min(abs(NASA_lat - .x))
    ),
    lon_idx = purrr::map_int(
      Longitude,
      ~ which.min(abs(NASA_lon - .x))
    )
  )


# ------------------------------------------------------------
# Read each NASA month ONCE
# ------------------------------------------------------------

NASA_results_list <- vector("list", length(NASA_files))

for (j in seq_along(NASA_files)) {
  
  cat(
    "Processing",
    j, "of", length(NASA_files),
    "-",
    format(nasa_years[j], "%Y-%m"),
    "\n"
  )
  
  h5f <- hdf5r::H5File$new(
    NASA_files[j],
    mode = "r"
  )
  
  precip <- h5f[["/Grid/precipitation"]]$read()
  
  h5f$close_all()
  
  # Extract precipitation for all islands
  values <- mapply(
    function(lat_idx, lon_idx) {
      value <- precip[lat_idx, lon_idx]
      
      if (value <= -9999) {
        value <- NA_real_
      }
      
      as.numeric(value)
    },
    NASA_indices$lat_idx,
    NASA_indices$lon_idx
  )
  
  # Convert mm/hour -> mm/month
  days <- lubridate::days_in_month(nasa_years[j])
  
  monthly_values <- values * 24 * days
  
  NASA_results_list[[j]] <- data.frame(
    date = nasa_years[j],
    year = lubridate::year(nasa_years[j]),
    month = lubridate::month(nasa_years[j]),
    monthly_precip = monthly_values,
    Island = NASA_indices$Island,
    product = "NASA"
  )
}


# ------------------------------------------------------------
# Combine all months
# ------------------------------------------------------------

NASA_monthly <- dplyr::bind_rows(NASA_results_list)


# ------------------------------------------------------------
# Diagnostic checks
# ------------------------------------------------------------

NASA_monthly %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_islands = dplyr::n_distinct(Island),
    first_date = min(date),
    last_date = max(date),
    n_missing = sum(is.na(monthly_precip)),
    min_monthly = min(monthly_precip, na.rm = TRUE),
    mean_monthly = mean(monthly_precip, na.rm = TRUE),
    max_monthly = max(monthly_precip, na.rm = TRUE)
  )

# ============================================================
# 9. ERA5
# ============================================================

# ============================================================
# ERA5 MONTHLY PRECIPITATION
# ============================================================

ERA5_path <- "C:/Users/Owner/OneDrive/Desktop/BAILEY/Precipitation/ERA5/5902d1318f4ade602cd74e19890b3ec1.nc"


# ------------------------------------------------------------
# 1. Check the ERA5 file
# ------------------------------------------------------------

nc <- ncdf4::nc_open(ERA5_path)

lat <- ncdf4::ncvar_get(nc, "latitude")
lon <- ncdf4::ncvar_get(nc, "longitude")

time_raw <- ncdf4::ncvar_get(nc, "valid_time")

time_units <- nc$dim[["valid_time"]]$units

tp_units <- nc$var[["tp"]]$units

ncdf4::nc_close(nc)

cat("ERA5 latitude range:", min(lat), "to", max(lat), "\n")
cat("ERA5 longitude range:", min(lon), "to", max(lon), "\n")
cat("Number of time steps:", length(time_raw), "\n")
cat("Time units:", time_units, "\n")
cat("TP units:", tp_units, "\n")


# ------------------------------------------------------------
# 2. Convert ERA5 time to dates
# ------------------------------------------------------------

ERA5_dates <- as.Date(
  as.POSIXct(
    time_raw,
    origin = "1970-01-01",
    tz = "UTC"
  )
)

cat("ERA5 date range:",
    format(min(ERA5_dates)),
    "to",
    format(max(ERA5_dates)),
    "\n")


# ------------------------------------------------------------
# 3. Keep only 1990-2024
# ------------------------------------------------------------

ERA5_keep <- ERA5_dates >= as.Date("1990-01-01") &
  ERA5_dates <= as.Date("2024-12-01")

ERA5_dates_use <- ERA5_dates[ERA5_keep]

ERA5_time_indices <- which(ERA5_keep)

cat("Number of ERA5 months used:",
    length(ERA5_dates_use),
    "\n")

cat("Analysis period:",
    format(min(ERA5_dates_use)),
    "to",
    format(max(ERA5_dates_use)),
    "\n")


# ------------------------------------------------------------
# 4. Function to extract ERA5 precipitation
# ------------------------------------------------------------

extract_and_aggregate_ERA5 <- function(
    target_lat,
    target_lon,
    island_name
) {
  
  # Open ERA5 file
  nc <- ncdf4::nc_open(ERA5_path)
  on.exit(
    ncdf4::nc_close(nc),
    add = TRUE
  )
  
  # Read coordinates
  lat <- ncdf4::ncvar_get(nc, "latitude")
  lon <- ncdf4::ncvar_get(nc, "longitude")
  
  # Find nearest grid cell
  lat_idx <- which.min(
    abs(lat - target_lat)
  )
  
  lon_idx <- which.min(
    abs(lon - target_lon)
  )
  
  # Extract ALL 1990-2024 monthly values
  #
  # tp dimensions are:
  # longitude x latitude x valid_time
  
  vec_raw <- ncdf4::ncvar_get(
    nc,
    "tp",
    start = c(
      lon_idx,
      lat_idx,
      min(ERA5_time_indices)
    ),
    count = c(
      1,
      1,
      length(ERA5_time_indices)
    )
  )
  
  days_per_month <- lubridate::days_in_month(ERA5_dates_use)
  
  monthly_precip_mm <-
    as.numeric(vec_raw) * 1000 * days_per_month
  
  
  # Build output data frame
  df <- data.frame(
    date = ERA5_dates_use,
    year = lubridate::year(ERA5_dates_use),
    month = lubridate::month(ERA5_dates_use),
    monthly_precip = monthly_precip_mm,
    Island = island_name,
    product = "ERA5"
  )
  
  
  # Check for impossible values
  df$monthly_precip[
    !is.finite(df$monthly_precip)
  ] <- NA_real_
  
  
  return(df)
}


# ------------------------------------------------------------
# 5. Test extraction for Laura before running all islands
# ------------------------------------------------------------

ERA5_test <- extract_and_aggregate_ERA5(
  target_lat = 7.15,
  target_lon = 171.03,
  island_name = "Laura"
)

cat(
  "Laura ERA5 rows:",
  nrow(ERA5_test),
  "\n"
)

cat(
  "Laura date range:",
  format(min(ERA5_test$date)),
  "to",
  format(max(ERA5_test$date)),
  "\n"
)

cat(
  "Laura missing months:",
  sum(is.na(ERA5_test$monthly_precip)),
  "\n"
)

print(head(ERA5_test))
print(tail(ERA5_test))


# ------------------------------------------------------------
# 6. Run extraction for all islands
# ------------------------------------------------------------

ERA5_results <- island_features %>%
  dplyr::mutate(
    ERA5 = purrr::pmap(
      list(
        Latitude,
        Longitude,
        Island
      ),
      ~ extract_and_aggregate_ERA5(
        target_lat = ..1,
        target_lon = ..2,
        island_name = ..3
      )
    )
  )


# ------------------------------------------------------------
# 7. Combine all islands
# ------------------------------------------------------------

ERA5_monthly <- dplyr::bind_rows(
  ERA5_results$ERA5
)


# ------------------------------------------------------------
# 8. Basic checks
# ------------------------------------------------------------

cat(
  "Total ERA5 rows:",
  nrow(ERA5_monthly),
  "\n"
)

cat(
  "Number of islands:",
  dplyr::n_distinct(ERA5_monthly$Island),
  "\n"
)

cat(
  "Date range:",
  format(min(ERA5_monthly$date)),
  "to",
  format(max(ERA5_monthly$date)),
  "\n"
)

cat(
  "Missing monthly values:",
  sum(is.na(ERA5_monthly$monthly_precip)),
  "\n"
)


# ------------------------------------------------------------
# 9. Check rows per island
# ------------------------------------------------------------

ERA5_monthly %>%
  dplyr::count(Island) %>%
  print(n = Inf)


# ------------------------------------------------------------
# 10. Check precipitation values
# ------------------------------------------------------------

summary(ERA5_monthly$monthly_precip)


# ------------------------------------------------------------
# 11. Save ERA5 results
# ------------------------------------------------------------

saveRDS(
  ERA5_monthly,
  "C:/Users/Owner/OneDrive/Desktop/BAILEY/Precipitation/ERA5/ERA5_monthly.rds"
)

# ============================================================
# 10. JAXA
# ============================================================


JAXA_dir <- file.path(precip_dir, "JAXA")

# ------------------------------------------------------------
# Find JAXA HDF5 files
# ------------------------------------------------------------

JAXA_files <- list.files(
  JAXA_dir,
  pattern = "\\.h5$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(JAXA_files) == 0) {
  stop("No JAXA HDF5 files found in: ", JAXA_dir)
}

cat("Number of JAXA files:", length(JAXA_files), "\n")
cat("First JAXA file:", basename(JAXA_files[1]), "\n")


# ------------------------------------------------------------
# Helper functions for JAXA grid
# ------------------------------------------------------------

nearest_lat_index <- function(grid_lat, target_lat) {
  
  if (!is.finite(target_lat)) {
    stop("Target latitude is not finite.")
  }
  
  which.min(abs(grid_lat - target_lat))
}


nearest_lon_index <- function(grid_lon, target_lon) {
  
  if (!is.finite(target_lon)) {
    stop("Target longitude is not finite.")
  }
  
  # JAXA longitude grid is -180 to 180, so make sure the
  # target longitude is represented in the same convention.
  target_lon <- ((target_lon + 180) %% 360) - 180
  
  which.min(abs(grid_lon - target_lon))
}


# ------------------------------------------------------------
# Read JAXA spatial grid
# ------------------------------------------------------------

h5_first <- hdf5r::H5File$new(
  JAXA_files[1],
  mode = "r"
)

JAXA_lat_raw <- h5_first[["/Grid/Latitude"]]$read()
JAXA_lon_raw <- h5_first[["/Grid/Longitude"]]$read()

h5_first$close_all()


# ------------------------------------------------------------
# Reconstruct 1-D coordinate axes
#
# IMPORTANT:
# JAXA Latitude is stored as a flattened 1,800 x 3,600
# spatial grid. Therefore length(JAXA_lat_raw) = 6,480,000.
#
# The latitude value is repeated across each longitude column.
# The longitude values represent the 3,600 longitude positions.
#
# We reconstruct the coordinate axes using unique values.
# ------------------------------------------------------------

JAXA_lat <- unique(as.numeric(JAXA_lat_raw))
JAXA_lon <- unique(as.numeric(JAXA_lon_raw))


# ------------------------------------------------------------
# Check reconstructed grid
# ------------------------------------------------------------

cat("\nJAXA grid check:\n")
cat("Raw latitude values:", length(JAXA_lat_raw), "\n")
cat("Raw longitude values:", length(JAXA_lon_raw), "\n")
cat("Unique latitude values:", length(JAXA_lat), "\n")
cat("Unique longitude values:", length(JAXA_lon), "\n")

cat(
  "Latitude range:",
  min(JAXA_lat, na.rm = TRUE),
  "to",
  max(JAXA_lat, na.rm = TRUE),
  "\n"
)

cat(
  "Longitude range:",
  min(JAXA_lon, na.rm = TRUE),
  "to",
  max(JAXA_lon, na.rm = TRUE),
  "\n"
)


# JAXA should contain a 0.1-degree global grid:
# 1,800 latitude positions
# 3,600 longitude positions

if (length(JAXA_lat) != 1800) {
  stop(
    "Unexpected JAXA latitude dimension: ",
    length(JAXA_lat),
    ". Expected 1800."
  )
}

if (length(JAXA_lon) != 3600) {
  stop(
    "Unexpected JAXA longitude dimension: ",
    length(JAXA_lon),
    ". Expected 3600."
  )
}


# ------------------------------------------------------------
# Check that the coordinates are ordered as expected
# ------------------------------------------------------------

if (any(diff(JAXA_lat) <= 0)) {
  stop("JAXA latitude coordinates are not strictly increasing.")
}

if (any(diff(JAXA_lon) <= 0)) {
  stop("JAXA longitude coordinates are not strictly increasing.")
}

cat(
  "JAXA latitude spacing:",
  unique(round(diff(JAXA_lat), 5))[1],
  "degrees\n"
)

cat(
  "JAXA longitude spacing:",
  unique(round(diff(JAXA_lon), 5))[1],
  "degrees\n"
)


# ------------------------------------------------------------
# Identify the JAXA grid cell corresponding to each island
# ------------------------------------------------------------

JAXA_indices <- island_features |>
  dplyr::mutate(
    lat_idx = purrr::map_int(
      Latitude,
      ~ nearest_lat_index(JAXA_lat, .x)
    ),
    
    lon_idx = purrr::map_int(
      Longitude,
      ~ nearest_lon_index(JAXA_lon, .x)
    )
  )


# Optional diagnostic showing the selected grid locations

JAXA_indices <- JAXA_indices |>
  dplyr::mutate(
    JAXA_latitude = JAXA_lat[lat_idx],
    JAXA_longitude = JAXA_lon[lon_idx],
    latitude_difference = JAXA_latitude - Latitude,
    longitude_difference = JAXA_longitude - Longitude
  )

cat("\nExample JAXA grid locations:\n")

print(
  JAXA_indices |>
    dplyr::select(
      Island,
      Latitude,
      Longitude,
      JAXA_latitude,
      JAXA_longitude
    ) |>
    head(10)
)


# ------------------------------------------------------------
# Function to extract one JAXA file
# ------------------------------------------------------------

extract_JAXA_file <- function(file_path, island_index_table) {
  
  # ----------------------------------------------------------
  # Determine date from filename
  #
  # JAXA filenames are formatted YYYYMMDD.h5
  # e.g. 19980101.h5
  # ----------------------------------------------------------
  
  file_name <- basename(file_path)
  
  date_string <- sub(
    "\\.h5$",
    "",
    file_name,
    ignore.case = TRUE
  )
  
  if (!grepl("^[0-9]{8}$", date_string)) {
    warning(
      "Could not determine date from JAXA filename: ",
      file_name
    )
    return(NULL)
  }
  
  file_date <- as.Date(
    date_string,
    format = "%Y%m%d"
  )
  
  if (is.na(file_date)) {
    warning(
      "Invalid date in JAXA filename: ",
      file_name
    )
    return(NULL)
  }
  
  
  # ----------------------------------------------------------
  # Open HDF5 file
  # ----------------------------------------------------------
  
  h5_file <- hdf5r::H5File$new(
    file_path,
    mode = "r"
  )
  
  on.exit(
    h5_file$close_all(),
    add = TRUE
  )
  
  
  # ----------------------------------------------------------
  # Read monthly precipitation rate
  # ----------------------------------------------------------
  
  precip <- h5_file[["/Grid/monthlyPrecipRate"]]$read()
  
  
  # ----------------------------------------------------------
  # Confirm expected dimensions
  # ----------------------------------------------------------
  
  precip_dim <- dim(precip)
  
  if (length(precip_dim) != 2) {
    stop(
      "Unexpected JAXA precipitation dimensions in ",
      file_name,
      ": ",
      paste(precip_dim, collapse = " x "),
      ". Expected a 2-D 1800 x 3600 matrix."
    )
  }
  
  if (precip_dim[1] != 1800 || precip_dim[2] != 3600) {
    stop(
      "Unexpected JAXA precipitation dimensions in ",
      file_name,
      ": ",
      paste(precip_dim, collapse = " x "),
      ". Expected 1800 x 3600."
    )
  }
  
  
  # ----------------------------------------------------------
  # Extract precipitation at each island
  #
  # monthlyPrecipRate is a 2-D matrix:
  #   rows    = latitude
  #   columns = longitude
  #
  # Use mapply so each island receives its corresponding
  # latitude/longitude pair.
  # ----------------------------------------------------------
  
  monthly_values <- mapply(
    FUN = function(lat_idx, lon_idx) {
      precip[lat_idx, lon_idx]
    },
    island_index_table$lat_idx,
    island_index_table$lon_idx
  )
  
  
  # ----------------------------------------------------------
  # Convert precipitation rate to monthly precipitation
  #
  # JAXA monthlyPrecipRate is in mm/hour.
  #
  # mm/hour × 24 hours/day × days/month
  # = mm/month
  # ----------------------------------------------------------
  
  days_this_month <- lubridate::days_in_month(
    file_date
  )
  
  monthly_values_mm <- as.numeric(monthly_values) *
    24 *
    days_this_month
  
  
  # ----------------------------------------------------------
  # Return island-level monthly values
  # ----------------------------------------------------------
  
  data.frame(
    date = file_date,
    year = lubridate::year(file_date),
    month = lubridate::month(file_date),
    island = island_index_table$Island,
    monthly_precip = monthly_values_mm,
    product = "JAXA",
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# Extract all JAXA files
# ------------------------------------------------------------

cat("\nExtracting JAXA precipitation...\n")

JAXA_results <- vector(
  "list",
  length(JAXA_files)
)

for (i in seq_along(JAXA_files)) {
  
  if (i %% 12 == 0 || i == 1 || i == length(JAXA_files)) {
    cat(
      "Processing JAXA file",
      i,
      "of",
      length(JAXA_files),
      "\n"
    )
  }
  
  JAXA_results[[i]] <- tryCatch(
    
    extract_JAXA_file(
      JAXA_files[i],
      JAXA_indices
    ),
    
    error = function(e) {
      
      warning(
        "Error processing ",
        basename(JAXA_files[i]),
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
}


# ------------------------------------------------------------
# Combine JAXA results
# ------------------------------------------------------------

JAXA_results_df <- dplyr::bind_rows(
  JAXA_results
)


# ------------------------------------------------------------
# Restrict analysis period
#
# JAXA begins in 1998.
# The overall analysis period is 1990–2024.
# ------------------------------------------------------------

JAXA_results_df <- JAXA_results_df[
  JAXA_results_df$year >= 1998 &
    JAXA_results_df$year <= 2024,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# Convert monthly JAXA precipitation to annual totals
#
# A year is retained only when all 12 months are available.
# ------------------------------------------------------------

JAXA_annual <- JAXA_results_df |>
  dplyr::group_by(
    island,
    year
  ) |>
  dplyr::summarise(
    n_months = sum(!is.na(monthly_precip)),
    annual_precip = ifelse(
      n_months == 12,
      sum(monthly_precip, na.rm = TRUE),
      NA_real_
    ),
    .groups = "drop"
  )


# Keep only complete years

JAXA_annual <- JAXA_annual[
  JAXA_annual$n_months == 12,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# Final JAXA annual dataset
# ------------------------------------------------------------

JAXA_annual <- JAXA_annual |>
  dplyr::select(
    island,
    year,
    annual_precip
  )

cat(
  "\nComplete JAXA annual observations:",
  nrow(JAXA_annual),
  "\n"
)

print(
  head(JAXA_annual)
)

# ============================================================
# 11. MSWEP
# ============================================================


MSWEP_dir <- file.path(precip_dir, "MSWEP")

# Find all MSWEP NetCDF files
MSWEP_files <- list.files(
  MSWEP_dir,
  pattern = "\\.nc$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(MSWEP_files) == 0) {
  stop("No MSWEP NetCDF files found in: ", MSWEP_dir)
}

cat("Number of MSWEP files:", length(MSWEP_files), "\n")


# ------------------------------------------------------------
# Read MSWEP spatial grid
# ------------------------------------------------------------

nc_first <- ncdf4::nc_open(MSWEP_files[1])

MSWEP_lon <- ncdf4::ncvar_get(nc_first, "lon")
MSWEP_lat <- ncdf4::ncvar_get(nc_first, "lat")

ncdf4::nc_close(nc_first)

cat("MSWEP longitude values:", length(MSWEP_lon), "\n")
cat("MSWEP latitude values:", length(MSWEP_lat), "\n")

cat(
  "MSWEP longitude range:",
  min(MSWEP_lon, na.rm = TRUE),
  "to",
  max(MSWEP_lon, na.rm = TRUE),
  "\n"
)

cat(
  "MSWEP latitude range:",
  min(MSWEP_lat, na.rm = TRUE),
  "to",
  max(MSWEP_lat, na.rm = TRUE),
  "\n"
)

# Confirm expected global 0.1-degree grid
if (length(MSWEP_lon) != 3600) {
  stop(
    "Unexpected MSWEP longitude dimension: ",
    length(MSWEP_lon),
    ". Expected 3600."
  )
}

if (length(MSWEP_lat) != 1800) {
  stop(
    "Unexpected MSWEP latitude dimension: ",
    length(MSWEP_lat),
    ". Expected 1800."
  )
}


# ------------------------------------------------------------
# Helper functions for nearest grid cell
# ------------------------------------------------------------

nearest_lat_index_MSWEP <- function(grid_lat, target_lat) {
  
  if (!is.finite(target_lat)) {
    stop("Target latitude is not finite.")
  }
  
  which.min(abs(grid_lat - target_lat))
}


nearest_lon_index_MSWEP <- function(grid_lon, target_lon) {
  
  if (!is.finite(target_lon)) {
    stop("Target longitude is not finite.")
  }
  
  # MSWEP uses -180 to 180 longitude
  target_lon <- ((target_lon + 180) %% 360) - 180
  
  which.min(abs(grid_lon - target_lon))
}


# ------------------------------------------------------------
# Determine MSWEP grid cell for every island
# ------------------------------------------------------------

MSWEP_indices <- island_features |>
  dplyr::mutate(
    lat_idx = purrr::map_int(
      Latitude,
      ~ nearest_lat_index_MSWEP(MSWEP_lat, .x)
    ),
    lon_idx = purrr::map_int(
      Longitude,
      ~ nearest_lon_index_MSWEP(MSWEP_lon, .x)
    )
  ) |>
  dplyr::mutate(
    MSWEP_latitude = MSWEP_lat[lat_idx],
    MSWEP_longitude = MSWEP_lon[lon_idx],
    latitude_difference = MSWEP_latitude - Latitude,
    longitude_difference = MSWEP_longitude - Longitude
  )

cat("\nMSWEP island/grid matching:\n")
print(
  MSWEP_indices |>
    dplyr::select(
      Island,
      Latitude,
      MSWEP_latitude,
      latitude_difference,
      Longitude,
      MSWEP_longitude,
      longitude_difference
    )
)


# ------------------------------------------------------------
# Extract one MSWEP file
# ------------------------------------------------------------

extract_MSWEP_file <- function(
    file_path,
    island_index_table
) {
  
  file_name <- basename(file_path)
  
  # Expected filenames: YYYYMMDD.nc
  date_string <- sub(
    "\\.nc$",
    "",
    file_name,
    ignore.case = TRUE
  )
  
  if (!grepl("^[0-9]{8}$", date_string)) {
    
    warning(
      "Skipping MSWEP file with unexpected filename: ",
      file_name
    )
    
    return(NULL)
  }
  
  file_date <- as.Date(
    date_string,
    format = "%Y%m%d"
  )
  
  if (is.na(file_date)) {
    
    warning(
      "Could not parse date from MSWEP file: ",
      file_name
    )
    
    return(NULL)
  }
  
  
  # Open file
  nc <- ncdf4::nc_open(file_path)
  
  on.exit(
    ncdf4::nc_close(nc),
    add = TRUE
  )
  
  
  # Read precipitation
  #
  # NetCDF variable is lon x lat x time,
  # but time has length 1 and ncvar_get()
  # drops that singleton dimension.
  #
  # Therefore:
  # dim(precip) = 3600 x 1800
  
  precip <- ncdf4::ncvar_get(
    nc,
    "precipitation"
  )
  
  
  # Confirm structure
  precip_dim <- dim(precip)
  
  if (
    length(precip_dim) != 2 ||
    precip_dim[1] != 3600 ||
    precip_dim[2] != 1800
  ) {
    
    stop(
      "Unexpected MSWEP precipitation dimensions in ",
      file_name,
      ": ",
      paste(precip_dim, collapse = " x "),
      ". Expected 3600 x 1800."
    )
  }
  
  
  # Extract the nearest grid cell for each island
  monthly_values <- vapply(
    seq_len(nrow(island_index_table)),
    function(j) {
      
      value <- precip[
        island_index_table$lon_idx[j],
        island_index_table$lat_idx[j]
      ]
      
      as.numeric(value)
    },
    numeric(1)
  )
  
  
  # MSWEP precipitation is already mm month-1.
  #
  # -9999 is the missing-value code and should be
  # converted to NA.
  
  monthly_values[
    monthly_values == -9999
  ] <- NA_real_
  
  
  data.frame(
    date = file_date,
    year = lubridate::year(file_date),
    month = lubridate::month(file_date),
    island = island_index_table$Island,
    monthly_precip = monthly_values,
    product = "MSWEP",
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# Process all MSWEP files
# ------------------------------------------------------------

MSWEP_results <- vector(
  "list",
  length(MSWEP_files)
)

for (i in seq_along(MSWEP_files)) {
  
  if (
    i %% 12 == 0 ||
    i == 1 ||
    i == length(MSWEP_files)
  ) {
    
    cat(
      "MSWEP:",
      i,
      "/",
      length(MSWEP_files),
      basename(MSWEP_files[i]),
      "\n"
    )
  }
  
  
  MSWEP_results[[i]] <- tryCatch(
    
    extract_MSWEP_file(
      MSWEP_files[i],
      MSWEP_indices
    ),
    
    error = function(e) {
      
      warning(
        "Error processing MSWEP file ",
        basename(MSWEP_files[i]),
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
}


# Combine all monthly results
MSWEP_results_df <- dplyr::bind_rows(
  MSWEP_results
)


# ------------------------------------------------------------
# Restrict analysis period
# ------------------------------------------------------------

MSWEP_results_df <- MSWEP_results_df[
  MSWEP_results_df$year >= 1990 &
    MSWEP_results_df$year <= 2019,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# Diagnostics
# ------------------------------------------------------------

cat(
  "\nMSWEP monthly records:",
  nrow(MSWEP_results_df),
  "\n"
)

cat(
  "MSWEP missing precipitation values:",
  sum(is.na(MSWEP_results_df$monthly_precip)),
  "\n"
)

cat(
  "MSWEP precipitation range:",
  min(
    MSWEP_results_df$monthly_precip,
    na.rm = TRUE
  ),
  "to",
  max(
    MSWEP_results_df$monthly_precip,
    na.rm = TRUE
  ),
  "mm/month\n"
)


MSWEP_month_counts <- MSWEP_results_df |>
  dplyr::group_by(island) |>
  dplyr::summarise(
    n_months = sum(
      !is.na(monthly_precip)
    ),
    first_date = min(
      date,
      na.rm = TRUE
    ),
    last_date = max(
      date,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(MSWEP_month_counts)


# ------------------------------------------------------------
# Calculate annual precipitation
# ------------------------------------------------------------

MSWEP_annual <- MSWEP_results_df |>
  dplyr::group_by(
    island,
    year
  ) |>
  dplyr::summarise(
    
    n_months = sum(
      !is.na(monthly_precip)
    ),
    
    annual_precip = ifelse(
      n_months == 12,
      sum(
        monthly_precip,
        na.rm = TRUE
      ),
      NA_real_
    ),
    
    .groups = "drop"
  )


# Retain only complete calendar years
MSWEP_annual <- MSWEP_annual[
  MSWEP_annual$n_months == 12,
  ,
  drop = FALSE
]


# Keep final columns
MSWEP_annual <- MSWEP_annual |>
  dplyr::select(
    island,
    year,
    annual_precip
  )


cat(
  "\nComplete MSWEP island-years:",
  nrow(MSWEP_annual),
  "\n"
)

print(
  head(MSWEP_annual)
)

# ============================================================
# 12. CDR
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "PROCESSING CDR\n"
)

cat(
  "========================================\n"
)


CDR_dir <-
  file.path(
    precip_dir,
    "CDR"
  )

CDR_files <-
  find_nc_files(
    CDR_dir
  )

cat(
  "CDR files found:",
  length(
    CDR_files
  ),
  "\n"
)


if (
  length(
    CDR_files
  ) == 0
) {
  
  warning(
    "No CDR NetCDF files found."
  )
  
  CDR_monthly <-
    NULL
  
} else {
  
  
  # ----------------------------------------------------------
  # Use first CDR file
  #
  # Your local file contains:
  #
  # precip[lon, lat, datetime]
  #
  # datetime:
  # days since 1990-01-01
  #
  # 408 monthly timestamps:
  # Jan 1990 - Dec 2023
  # ----------------------------------------------------------
  
  CDR_file <-
    CDR_files[1]
  
  
  cat(
    "Using CDR file:",
    CDR_file,
    "\n"
  )
  
  
  nc <-
    ncdf4::nc_open(
      CDR_file
    )
  
  
  CDR_lat <-
    ncdf4::ncvar_get(
      nc,
      "lat"
    )
  
  CDR_lon <-
    ncdf4::ncvar_get(
      nc,
      "lon"
    )
  
  CDR_time <-
    ncdf4::ncvar_get(
      nc,
      "datetime"
    )
  
  
  CDR_dates <-
    as.Date(
      CDR_time,
      origin =
        "1990-01-01"
    )
  
  
  # ----------------------------------------------------------
  # Verify expected date range
  # ----------------------------------------------------------
  
  expected_CDR_dates <-
    seq.Date(
      from =
        as.Date(
          "1990-01-01"
        ),
      to =
        as.Date(
          "2023-12-01"
        ),
      by =
        "month"
    )
  
  
  if (
    length(CDR_dates) != length(expected_CDR_dates) ||
    any(CDR_dates != expected_CDR_dates)
  ) {
    warning(
      "CDR dates do not match expected Jan 1990-Dec 2023 sequence."
    )
  }
  
  
  # ----------------------------------------------------------
  # Find precipitation variable
  # ----------------------------------------------------------
  
  CDR_precip <-
    ncdf4::ncvar_get(
      nc,
      "precip"
    )
  
  
  ncdf4::nc_close(
    nc
  )
  
  
  # ----------------------------------------------------------
  # Find grid indices
  # ----------------------------------------------------------
  
  CDR_indices <-
    island_features %>%
    
    dplyr::mutate(
      
      lat_idx =
        purrr::map_int(
          Latitude,
          ~ nearest_lat_index(
            CDR_lat,
            .x
          )
        ),
      
      lon_idx =
        purrr::map_int(
          Longitude,
          ~ nearest_lon_index(
            CDR_lon,
            .x
          )
        )
    )
  
  
  # ----------------------------------------------------------
  # Extract CDR
  # ----------------------------------------------------------
  
  CDR_list <-
    vector(
      "list",
      nrow(
        CDR_indices
      )
    )
  
  
  for (
    j in seq_len(
      nrow(
        CDR_indices
      )
    )
  ) {
    
    island_values <-
      CDR_precip[
        CDR_indices$lon_idx[j],
        CDR_indices$lat_idx[j],
        ,
        drop = TRUE
      ]
    
    
    island_values <-
      as.numeric(
        island_values
      )
    
    
    # --------------------------------------------------------
    # CDR missing-value code:
    # -99
    # --------------------------------------------------------
    
    island_values[
      island_values == -99
    ] <-
      NA_real_
    
    
    CDR_list[[j]] <-
      data.frame(
        
        Island =
          CDR_indices$Island[j],
        
        Latitude =
          CDR_indices$Latitude[j],
        
        Longitude =
          CDR_indices$Longitude[j],
        
        date =
          CDR_dates,
        
        monthly_precip =
          island_values
      )
  }
  
  
  # ----------------------------------------------------------
  # CDR file is already monthly.
  # No conversion.
  # ----------------------------------------------------------
  
  CDR_monthly <-
    dplyr::bind_rows(
      CDR_list
    ) %>%
    
    dplyr::mutate(
      
      year =
        lubridate::year(
          date
        ),
      
      month =
        lubridate::month(
          date
        ),
      
      product =
        "CDR"
    ) %>%
    
    dplyr::select(
      Island,
      Latitude,
      Longitude,
      date,
      year,
      month,
      monthly_precip,
      product
    )
  
  
  saveRDS(
    CDR_monthly,
    file.path(
      rds_dir,
      "CDR_monthly.rds"
    )
  )
}


# ============================================================
# 13. COMBINE ALL PRODUCTS
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "COMBINING PRODUCTS\n"
)

cat(
  "========================================\n"
)


JAXA_monthly <- JAXA_results_df
MSWEP_monthly <- MSWEP_results_df


JAXA_monthly <- JAXA_results_df |>
  dplyr::rename(
    Island = island
  ) |>
  dplyr::left_join(
    island_features |>
      dplyr::select(
        Island,
        Latitude,
        Longitude
      ),
    by = "Island"
  )

MSWEP_monthly <- MSWEP_results_df |>
  dplyr::rename(
    Island = island
  ) |>
  dplyr::left_join(
    island_features |>
      dplyr::select(
        Island,
        Latitude,
        Longitude
      ),
    by = "Island"
  )



precip_monthly <-
  dplyr::bind_rows(
    NASA_monthly,
    ERA5_monthly,
    JAXA_monthly,
    MSWEP_monthly,
    CDR_monthly
  )


# ============================================================
# 14. KEEP REQUESTED DATE RANGES
# ============================================================

precip_monthly |>
  dplyr::group_by(product) |>
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_islands = dplyr::n_distinct(Island),
    n_missing = sum(is.na(monthly_precip)),
    first_date = min(date, na.rm = TRUE),
    last_date = max(date, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 15. CHECK FOR DUPLICATES
# ============================================================

duplicates <-
  precip_monthly %>%
  
  dplyr::count(
    Island,
    product,
    date
  ) %>%
  
  dplyr::filter(
    n > 1
  )


if (
  nrow(
    duplicates
  ) > 0
) {
  
  warning(
    "Duplicate island/product/date observations detected."
  )
  
  print(
    duplicates
  )
}


# ============================================================
# 16. MONTHLY SUMMARY
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "MONTHLY SUMMARY\n"
)

cat(
  "========================================\n"
  
  
)


monthly_summary <-
  precip_monthly %>%
  
  dplyr::group_by(
    product
  ) %>%
  
  dplyr::summarise(
    
    n_rows =
      dplyr::n(),
    
    n_islands =
      dplyr::n_distinct(
        Island
      ),
    
    first_date =
      min(
        date,
        na.rm = TRUE
      ),
    
    last_date =
      max(
        date,
        na.rm = TRUE
      ),
    
    n_missing =
      sum(
        is.na(
          monthly_precip
        )
      ),
    
    mean_monthly_precip =
      mean(
        monthly_precip,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  )


print(
  monthly_summary
)


# ============================================================
# 17. CALCULATE ANNUAL PRECIPITATION
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "CALCULATING ANNUAL PRECIPITATION\n"
)

cat(
  "========================================\n"
)


annual_precip <-
  precip_monthly %>%
  
  dplyr::group_by(
    
    Island,
    
    Latitude,
    
    Longitude,
    
    product,
    
    year
  ) %>%
  
  dplyr::summarise(
    
    n_months =
      dplyr::n_distinct(
        month[
          !is.na(
            monthly_precip
          )
        ]
      ),
    
    annual_precip_mm =
      
      ifelse(
        
        n_months == 12,
        
        sum(
          monthly_precip,
          na.rm = TRUE
        ),
        
        NA_real_
      ),
    
    .groups =
      "drop"
  )


# ============================================================
# 18. SAVE MONTHLY DATA
# ============================================================

monthly_csv <-
  file.path(
    output_dir,
    "precipitation_monthly_long.csv"
  )

monthly_rds <-
  file.path(
    rds_dir,
    "precipitation_monthly_long.rds"
  )


readr::write_csv(
  precip_monthly,
  monthly_csv
)


saveRDS(
  precip_monthly,
  monthly_rds
)


# ============================================================
# 19. SAVE ANNUAL DATA
# ============================================================

annual_csv <-
  file.path(
    output_dir,
    "precipitation_annual.csv"
  )

annual_rds <-
  file.path(
    rds_dir,
    "precipitation_annual.rds"
  )


readr::write_csv(
  annual_precip,
  annual_csv
)


saveRDS(
  annual_precip,
  annual_rds
)


# ============================================================
# 20. FINAL SUMMARY
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "SCRIPT COMPLETE\n"
)

cat(
  "========================================\n"
)


cat(
  "\nProducts processed:\n"
)

print(
  unique(
    precip_monthly$product
  )
)


cat(
  "\nNumber of islands:",
  dplyr::n_distinct(
    precip_monthly$Island
  ),
  "\n"
)


cat(
  "\nNumber of monthly observations:",
  nrow(
    precip_monthly
  ),
  "\n"
)


cat(
  "\nNumber of complete annual observations:",
  sum(
    !is.na(
      annual_precip$annual_precip_mm
    )
  ),
  "\n"
)


cat(
  "\nMonthly CSV:\n",
  monthly_csv,
  "\n"
)


cat(
  "\nAnnual CSV:\n",
  annual_csv,
  "\n"
)

