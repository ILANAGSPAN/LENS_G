#### More figures for the appendix?

library(tidyr)
library(dplyr)
library(patchwork)
library(rnaturalearth)
library(ggplot2)
library(lubridate)
library(purrr)
library(dplyr)
library(grid)

setwd("Q:/Hydro/projects/indopacom_PRIVATE/ILANA")  # off-site path directory

################################################################################
#PLOTTING-----------------------------------------------------------------------
################################################################################


#MAP OF REMOTE SENSING PRODUCT GRIDS
world <- ne_countries(scale = "medium", returnclass = "sf")

precip_grid <- ggplot(data = world) +
  geom_sf(fill = "lightgray", color = "black") +
  annotate("rect", xmin = 114.125, xmax = 114.375, ymin = 10.875, ymax = 11.125,
           color = "red", fill = NA, size = 1.2) +
  annotate("rect", xmin = 113.875, xmax = 114.125, ymin = 10.875, ymax = 11.125,
           color = "red", fill = NA, size = 1.2) +
  annotate("rect", xmin = 112.875, xmax = 113.125, ymin = 9.375, ymax = 9.625,
           color = "red", fill = NA, size = 1.2) +
  annotate("rect", xmin = 114.2, xmax = 114.3, ymin = 11, ymax = 11.1,
           color = "purple", fill = NA, size = 1.2) +
  annotate("rect", xmin = 114, xmax = 114.1, ymin = 10.9, ymax = 11,
           color = "purple", fill = NA, size = 1.2) +
  annotate("rect", xmin = 112.8, xmax = 112.9, ymin = 9.5, ymax = 9.6,
           color = "purple", fill = NA, size = 1.2) +
  coord_sf(xlim = c(112, 115), ylim = c(9, 12)) +
  geom_point(aes(x = 114.283, y = 11.05), color = "black", size = 1.5) +
  geom_point(aes(x = 114.085, y = 10.92), color = "black", size = 1.5) +
  geom_point(aes(x = 112.89, y = 9.55), color = "black", size = 1.5) +
  labs(x = "Longitude",
       y = "Latitude") +
  theme_minimal()

ggsave("R_Plots/Precipitation_Grid.png",
       precip_grid,
       width  = 15,
       height = 15,   # nudged down from 16 — the shorter subtitle frees vertical space
       units  = "cm",
       dpi    = 600)


ggsave("R_Plots/Precipitation_Grid.pdf",
       precip_grid,
       width  = 15,
       height = 15,   # nudged down from 16 — the shorter subtitle frees vertical space
       units  = "cm",
       dpi    = 600)

#-------------------------------------------------------------------------------


kwajalein <- joint_precip_df %>%
  filter(island == "Kwajelein")

kwajalein$date <- as.Date(kwajalein$date)

kwajalein_2010_2019 <- kwajalein %>%
  filter(date >= as.Date("2010-01-01"),
         date <= as.Date("2019-12-31"))


kwajalein_precip <- ggplot(
  kwajalein_2010_2019,
  aes(x = date,
      y = monthly_precip,
      color = product)
) +
  
  geom_line(
    linewidth = 0.5,
    alpha = 0.9
  ) +
  
  scale_color_manual(values = c(
    "CDR"   = "#D55E00",   # Vermillion
    "ERA5"  = "#009E73",   # Bluish Green
    "JAXA"  = "#0072B2",   # Blue
    "MSWEP" = "#E69F00",   # Orange
    "NASA"  = "#CC79A7"    # Reddish Purple
  )) +
  
#  scale_x_date(
#    date_breaks = "1 year",
#    date_labels = "%Y",
 #   expand = expansion(mult = c(0.01, 0.01))
#  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.03))
  ) +
  
  labs(
    x = "Year",
    y = "Monthly precipitation (mm)",
    color = NULL
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.justification = "top",
    
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    
    legend.key.size = unit(0.3, "cm"),
    legend.spacing = unit(0.05, "cm"),   # single spacing argument
    
    legend.box.margin = margin(0, 0, 0, -10),
    legend.margin = margin(0, 0, 0, 0),
    
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    
    plot.margin = margin(2, 2, 2, 2)
  )

ggsave("R_Plots/kwaj_Precip.png",
       kwajalein_precip,
       width  = 15,
       height = 7,   # nudged down from 16 — the shorter subtitle frees vertical space
       units  = "cm",
       dpi    = 600)


ggsave("R_Plots/kwaj_Precip.pdf",
       kwajalein_precip,
       width  = 15,
       height = 7,   # nudged down from 16 — the shorter subtitle frees vertical space
       units  = "cm",
       dpi    = 600)


#AVERAGE MONTHLY PRECIPITATION--------------------------------------------------

kwajalein_2000_2019 <- kwajalein %>%
  filter(date >= as.Date("2000-01-01"),
         date <= as.Date("2019-12-31"))

avg_monthly <- kwajalein_2000_2019 %>%
  filter(monthly_precip > 0) %>%
  mutate(month = month(date)) %>%
  group_by(product, month) %>%
  summarize(avg_precip = mean(monthly_precip)) %>%
  arrange(product, month)

# Convert month number to abbreviation
avg_monthly$month <- factor(
  avg_monthly$month,
  levels = 1:12,
  labels = month.abb,
  ordered = TRUE
)


avg_monthly_plot <- ggplot(
  avg_monthly,
  aes(x = month,
      y = avg_precip,
      color = product,
      group = product)
) +
  
  geom_line(
    linewidth = 0.5,
    alpha = 0.9
  ) +
  
  scale_color_manual(values = c(
    "CDR"   = "#D55E00",   # Vermillion
    "ERA5"  = "#009E73",   # Bluish Green
    "JAXA"  = "#0072B2",   # Blue
    "MSWEP" = "#E69F00",   # Orange
    "NASA"  = "#CC79A7"    # Reddish Purple
  )) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.03))
  ) +
  
  labs(
    x = "Month",
    y = "Average precipitation (mm)",
    color = NULL
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.justification = "top",
    
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    
    legend.key.size = unit(0.3, "cm"),
    legend.spacing = unit(0.05, "cm"),   # single spacing argument
    
    legend.box.margin = margin(0, 0, 0, -10),
    legend.margin = margin(0, 0, 0, 0),
    
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    
    plot.margin = margin(2, 2, 2, 2)
  )

ggsave("R_Plots/avg_monthly_plot.png",
       avg_monthly_plot,
       width  = 15,
       height = 7,   # nudged down from 16 — the shorter subtitle frees vertical space
       units  = "cm",
       dpi    = 600)


ggsave("R_Plots/avg_monthly_plot.pdf",
       avg_monthly_plot,
       width  = 15,
       height = 7,   # nudged down from 16 — the shorter subtitle frees vertical space
       units  = "cm",
       dpi    = 600)



################################################################################
#PLOTTING RECHARGE AND WIDTH BY REGION AND PRODUCT------------------------------
################################################################################


island_features <- read.csv("~/Spratly Features.csv")


islands_recharge_long <- island_features %>%
  filter(Island.ID <= 41 | Island.ID >= 48) %>%
  pivot_longer(cols = MSWEP:ERA5,
               names_to = "Product",
               values_to = "Recharge")


ggplot(islands_recharge_long, aes(x = Recharge, color = Product, fill = Product)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~Region) +
  ylab("Density")

ggsave(filename = "R_Plots/Recharge_by_Region.png", width = 10, height = 3, dpi = 300)


#width-------------------------------------------------

island_width_long <- pivot_longer(islands_recharge_long,
                                  cols = Width.1:Width.4,
                                  names_to = "Trial",
                                  values_to = "Width") 

ggplot(island_width_long, aes(x = Width)) +
  geom_histogram() +
  facet_wrap(~Region)


ggplot(island_width_long, aes(x = Width)) +
  geom_density() +
  xlab("Island Width (m)") +
  ylab("Density") +
  facet_wrap(~Region)

ggsave(filename = "R_Plots/Width_by_Region.png", width = 10, height = 3, dpi = 300)

#-------------------------------------------------------------------------------

library(tidyverse)

# Merge everything and filter for Laura
laura_precip_df <- bind_rows(
  list_rbind(CDR_results$CDR),
  list_rbind(NASA_results$NASA),
  list_rbind(JAXA_results$JAXA),
  list_rbind(MSWEP_results$MSWEP),
  list_rbind(ERA5_results$ERA5)
) %>%
  rename(island = island) %>% 
  # Filter down to only keep Laura island data
  filter(island == "Laura") %>% 
  dplyr::select(date, island, monthly_precip, product)



ggplot(laura_precip_df, aes(x = date, y = monthly_precip, color = product)) +
  # Draw the time-series lines
  geom_line(alpha = 0.8, size = 0.8) + 
  
  # Clean, publication-ready design theme
  theme_minimal() +
  labs(
    title = "Monthly Precipitation Timeline for Laura Island",
    subtitle = "Comparison across satellite and reanalysis data products",
    x = "Date",
    y = "Precipitation (mm / month)",
    color = "Data Product"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )



#RECREATE THIS BY FEEDING AVERAGE FOR PRECIPITATION AND WIDTH AND VARYING HYDRAULIC CAPCITY
#NASA
#---------------------------------------------------------------------------------------

library(hdf5r)
library(dplyr)
library(reshape2)
library(ggplot2)

path <- "Q:/Hydro/projects/indopacom_PRIVATE/ILANA/Data_Files/Precipitation/NASA"
files <- list.files(path, pattern = "\\.HDF5$", full.names = TRUE)

# Open first file to get lat/lon dimensions
h5 <- H5File$new(files[1], mode = "r")
lat <- h5[["/Grid/lat"]]$read()
lon <- h5[["/Grid/lon"]]$read()
nlat <- length(lat)
nlon <- length(lon)
h5$close_all()

# Initialize sum and count matrices
precip_sum <- matrix(0, nrow = nlat, ncol = nlon)
precip_count <- matrix(0, nrow = nlat, ncol = nlon)

# Loop through files, read one at a time
for (f in files) {
  cat("Processing:", f, "\n")
  h5 <- H5File$new(f, mode = "r")
  precip <- h5[["/Grid/precipitation"]]$read()
  h5$close_all()
  # Filter out invalid data
  precip[precip < 0] <- NA
  # Add valid data to running sum and count
  valid <- !is.na(precip)
  precip_sum[valid] <- precip_sum[valid] + precip[valid]
  precip_count[valid] <- precip_count[valid] + 1
}

# Compute mean (avoiding division by zero)
mean_precip <- precip_sum / precip_count
mean_precip[is.na(mean_precip)] <- 0


precip_df <- melt(mean_precip)
names(precip_df) <- c("lat_idx", "lon_idx", "mean_precip")
precip_df$lat <- lat[precip_df$lat_idx]
precip_df$lon <- lon[precip_df$lon_idx]

lat_min <- 3
lat_max <- 18
lon_min <- 100
lon_max <- 170

region_df <- precip_df %>%
  filter(lat >= lat_min, lat <= lat_max, lon >= lon_min, lon <= lon_max)


#JAXA=============================================================


# --- File path ---
path_JAXA <- "Q:/Hydro/projects/indopacom_PRIVATE/ILANA/Data_Files/Precipitation/JAXA"
JAXA_files <- list.files(path_JAXA, pattern = "\\.h5$", full.names = TRUE)

if (length(JAXA_files) == 0) stop("No JAXA .h5 files found.")

# --- Open first file to detect structure ---
h5 <- H5File$new(JAXA_files[1], mode = "r")

# Read lat/lon (2D arrays)
lat_mat <- h5[["/Grid/Latitude"]]$read()
lon_mat <- h5[["/Grid/Longitude"]]$read()

# Convert to 1D vectors
lat <- lat_mat[,1]
lon <- lon_mat[1,]

nlat <- length(lat)
nlon <- length(lon)

# --- Detect precipitation dataset inside /Grid ---
grid_items <- h5[["Grid"]]$ls()$name

candidate_datasets <- c(
  "monthlyPrecipRateGC",   # gauge-corrected monthly
  "monthlyPrecipRate",     # uncorrected monthly
  "precipitation",
  "precipRate"
)

precip_name <- intersect(candidate_datasets, grid_items)[1]

if (is.na(precip_name)) {
  stop("No known precipitation dataset found in /Grid.")
}

cat("Using precipitation dataset:", precip_name, "\n")

# --- Direct attribute access (Windows-safe) ---
dset <- h5[[paste0("/Grid/", precip_name)]]

units <- NA
if (dset$attr_exists("Units")) {
  units <- dset$attr[["Units"]]
} else if (dset$attr_exists("units")) {
  units <- dset$attr[["units"]]
}

cat("Units detected:", units, "\n")

h5$close_all()

# --- Initialize accumulators ---
precip_sum <- matrix(0, nrow = nlat, ncol = nlon)
precip_count <- matrix(0, nrow = nlat, ncol = nlon)
# --- Loop through files ---
failed_files <- character(0)

for (f in JAXA_files) {
  
  cat("Processing:", basename(f), "\n")
  
  tryCatch({
    
    h5 <- H5File$new(f, mode = "r")
    
    precip <- h5[[paste0("/Grid/", precip_name)]]$read()
    
    h5$close_all()
    
    # Dimension check
    if (!all(dim(precip) == c(nlat, nlon))) {
      stop("Dimension mismatch")
    }
    
    # Invalid values
    precip[precip < 0] <- NA
    
    # Accumulate
    valid <- !is.na(precip)
    precip_sum[valid] <- precip_sum[valid] + precip[valid]
    precip_count[valid] <- precip_count[valid] + 1
    
  }, error = function(e) {
    
    cat("  --> Skipping", basename(f), "\n")
    cat("      Reason:", conditionMessage(e), "\n")
    
    failed_files <<- c(failed_files, basename(f))
    
    # Close file if it was opened
    try(h5$close_all(), silent = TRUE)
    
  })
}

# --- Compute long-term mean ---
mean_precip <- precip_sum / precip_count
mean_precip[is.na(mean_precip)] <- 0

# --- Convert to tidy data frame ---
precip_JAXA <- melt(mean_precip)
names(precip_JAXA) <- c("lat_idx", "lon_idx", "mean_precip")

precip_JAXA$lat <- lat[precip_JAXA$lat_idx]
precip_JAXA$lon <- lon[precip_JAXA$lon_idx]

# --- Subset region ---
lat_min <- 3
lat_max <- 18
lon_min <- 100
lon_max <- 170

region_JAXA <- precip_JAXA %>%
  filter(lat >= lat_min, lat <= lat_max,
         lon >= lon_min, lon <= lon_max)


#MSWEP===========================================================================

library(hdf5r)
library(reshape2)
library(dplyr)
library(ggplot2)
library(viridis)
library(lubridate)
library(terra)
library(dplyr)

# path to your MSWEP files
path_mswep <- "Q:/Hydro/projects/indopacom_PRIVATE/ILANA/Data_Files/Precipitation/MSWEP"

# list all NetCDF files
mswep_files <- list.files(path_mswep, pattern = "\\.nc$", full.names = TRUE)


# Create a SpatRaster from all files
r_mswep <- rast(mswep_files)

# Check variable names
names(r_mswep)

mean_precip_mswep <- mean(r_mswep, na.rm = TRUE)


lat_min <- 3
lat_max <- 18
lon_min <- 100
lon_max <- 170

region_mswep <- crop(mean_precip_mswep, ext(lon_min, lon_max, lat_min, lat_max))
plot(region_mswep)

precip_mswep <- as.data.frame(region_mswep, xy = TRUE)
names(precip_mswep) <- c("lon", "lat", "mean_precip")


#CDR----------------------------------------------------------------------------


path_CDR <- "Q:/Hydro/projects/indopacom_PRIVATE/ILANA/Data_Files/Precipitation/CDR/CDR_2025-07-15114010am"

file <- list.files(
  path_CDR,
  pattern = "\\.nc$",
  full.names = TRUE
)

nc <- nc_open(file)

lat <- ncvar_get(nc, "lat")
lon <- ncvar_get(nc, "lon")

# precip[lon, lat, time]
precip <- ncvar_get(nc, "precip")

# Replace fill values with NA
fill <- ncatt_get(nc, "precip", "_FillValue")$value
precip[precip == fill] <- NA

# Mean across all months
mean_precip <- apply(precip, c(1, 2), mean, na.rm = TRUE)

# Convert from [lon, lat] to [lat, lon]
mean_precip <- t(mean_precip)

nc_close(nc)

# Convert to data frame
precip_CDR <- melt(mean_precip)
names(precip_CDR) <- c("lat_idx", "lon_idx", "mean_precip")

precip_CDR$lat <- lat[precip_CDR$lat_idx]
precip_CDR$lon <- lon[precip_CDR$lon_idx]

# Crop study region
region_CDR <- precip_CDR %>%
  filter(
    lat >= 3,
    lat <= 18,
    lon >= 100,
    lon <= 170
  )


#ERA5----------------------------------------------------------------------------

library(terra)
library(dplyr)
library(ggplot2)

path_era5 <- "Q:/Hydro/projects/indopacom_PRIVATE/ILANA/Data_Files/Precipitation/ERA5"
era5_files <- list.files(path_era5, pattern = "\\.nc$", full.names = TRUE)


library(terra)


path_ERA5 <- "Q:/Hydro/projects/indopacom_PRIVATE/ILANA/Data_Files/Precipitation/ERA5/1b6b56dac5acfe0ed774944eadf37725/"

era_files <- list.files(
  path_ERA5,
  pattern = "\\.nc$",
  full.names = TRUE
)

# Read both NetCDF files
r <- rast(era_files)

# Keep only total precipitation layers
tp <- r[[grep("^tp_", names(r))]]

# Compute long-term mean
mean_precip_era5 <- mean(tp, na.rm = TRUE)

# Convert meters --> millimeters
mean_precip_era5 <- mean_precip_era5 * 1000

# Crop study area
region_era5 <- crop(
  mean_precip_era5,
  ext(100, 170, 3, 18)
)

# Convert to data frame
precip_ERA5 <- as.data.frame(region_era5, xy = TRUE, na.rm = TRUE)
names(precip_ERA5) <- c("lon", "lat", "mean_precip")



#==============================================================================
# NASA
#==============================================================================

library(ggplot2)
library(patchwork)
library(viridis)
library(grid)

#==============================================================================
# Common theme
#==============================================================================

base_theme <- theme_bw(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.4),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    plot.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.height = unit(0.7, "cm"),
    legend.key.width = unit(0.35, "cm")
  )

#==============================================================================
# Shared precipitation color scale
#==============================================================================

precip_scale <- scale_fill_viridis_c(
  option = "C",
  limits = c(0, 600),
  breaks = seq(0, 600, 200),
  oob = scales::squish,
  name = "Precipitation\n(mm month⁻¹)"
)

#==============================================================================
# Region boxes
#==============================================================================

add_region_boxes <- function(p) {
  
  p +
    # Spratly Islands
    annotate(
      "rect",
      xmin = 111.9,
      xmax = 114.0,
      ymin = 8.6,
      ymax = 10.0,
      fill = NA,
      colour = "black",
      linewidth = 0.45
    ) +
    annotate(
      "label",
      x = 112.1,
      y = 12.6,
      label = "Spratlys",
      hjust = 0,
      size = 2.5,
      fontface = "bold",
      fill = "white"
    ) +
    
    # Micronesia
    annotate(
      "rect",
      xmin = 137,
      xmax = 163,
      ymin = 4,
      ymax = 11,
      fill = NA,
      colour = "black",
      linewidth = 0.45
    ) +
    annotate(
      "label",
      x = 137.3,
      y = 12.6,
      label = "Micronesia",
      hjust = 0,
      size = 2.5,
      fontface = "bold",
      fill = "white"
    )
}

#==============================================================================
# Helper function for all precipitation maps
#==============================================================================

make_precip_plot <- function(data,
                             subtitle,
                             fill_var,
                             show_x = FALSE,
                             show_y = TRUE,
                             title = NULL){
  
  p <- ggplot(
    data,
    aes(
      x = lon,
      y = lat,
      fill = {{ fill_var }}
    )
  ) +
    
    geom_raster() +
    
    precip_scale +
    
    coord_equal() +
    
    labs(
      title = title,
      subtitle = subtitle,
      x = if(show_x) "Longitude" else NULL,
      y = if(show_y) "Latitude" else NULL
    ) +
    
    base_theme
  
  if(!show_x){
    p <- p +
      theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank()
      )
  }
  
  if(!show_y){
    p <- p +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()
      )
  }
  
  add_region_boxes(p)
}

#==============================================================================
# Individual panels
#==============================================================================

NASA <- make_precip_plot(
  data     = region_df,
  subtitle = "NASA GPM IMERG",
  fill_var = mean_precip * 744,
  show_x   = FALSE,
  show_y   = TRUE,
  title    = "30-Year Mean Monthly Precipitation"
)

JAXA <- make_precip_plot(
  data     = region_JAXA,
  subtitle = "JAXA GSMaP",
  fill_var = mean_precip * 744,
  show_x   = FALSE,
  show_y   = TRUE
)

MSWEP <- make_precip_plot(
  data     = precip_mswep,
  subtitle = "MSWEP",
  fill_var = mean_precip,
  show_x   = TRUE,     # now the bottom panel — needs the x-axis
  show_y   = TRUE
)

#==============================================================================
# Combine figure — single column, three stacked panels
#==============================================================================

precip_products <-
  NASA / JAXA / MSWEP +
  plot_layout(
    guides  = "collect",
    heights = c(1, 1, 1)
  ) +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "right",
    plot.tag        = element_text(face = "bold", size = 10)
  )

precip_products

#==============================================================================
# Save
#==============================================================================

ggsave(
  filename = "Q:/Hydro/projects/indopacom_PRIVATE/ILANA/R_Plots/Figure_Precipitation_Comparison.png",
  plot     = precip_products,
  width    = 15,
  height   = 13,   # shorter now — only 3 stacked panels instead of a 3x2 grid
  units    = "cm",
  dpi      = 600
)

ggsave(
  filename = "R_Plots/Figure_Precipitation_Comparison.pdf",
  plot     = precip_products,
  width    = 15,
  height   = 13,
  units    = "cm"
)
#-------------------------------------------------------------------------------


ggplot(final_results, aes(x = Island.ID, y = Z_MAX_avg_precip)) +
  geom_errorbar(aes(x = Island.ID - 0.2, ymin = Z_MAX_CDR, ymax = Z_MAX_NASA, color = "Precipitation"), width = 0.1, linewidth = 1) +
  geom_errorbar(aes(ymin = Z_MAX_width1, ymax = Z_MAX_width4, color = "Width"), width = 0.2, linewidth = 1) +
  geom_errorbar(aes(x = Island.ID + 0.2, ymin = Z_MAX_high_hc, ymax = Z_MAX_low_hc, color = "Hydraulic Conductivity"), width = 0.2, linewidth = 1) +
  geom_point(aes(shape = Region), size = 3.5, color = "black") +
  scale_color_manual(name = "Uncertainty source", values = c("Precipitation" = "blue", "Width" = "green", "Hydraulic Conductivity" = "red")) +
#  scale_shape_manual(name = "Region",
 #   values = c(
#      "FSM" = 16,   # filled circle
#      "RMI" = 17,   # triangle
 #     "Kiribati" = 15)) +
  labs(x = "Island ID",
       y = "Predicted Z Max")


ggsave(filename = "Ranges.png", width = 14, height = 6, dpi = 300)

#--------------------------------------------------------------------------------

#RAINFALL ANALYSIS OCTOBER 2016
NASA = ggplot(region_df, aes(x = lon, y = lat, fill = precip)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C",
                       limits = c(0, 1100),
                       breaks = c(0, 200, 400, 600, 800, 1000)) +
  coord_fixed() +
  labs(subtitle = "NASA GPM IMERG") +
  theme_minimal()

ERA5 = ggplot(region_era5, aes(x = x, y = y, fill = precipitation)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C",
                       limits = c(0, 1000),
                       breaks = c(0, 200, 400, 600, 800, 1000)) +
  coord_fixed() +
  labs(subtitle = "ERA5") +
  theme_minimal()

CDR = ggplot(region_cdr, aes(x = x, y = y, fill = precip)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C",
                       limits = c(0, 1000),
                       breaks = c(0, 200, 400, 600, 800, 1000)) +
  coord_fixed() +
  labs(subtitle = "PERSIANN-CDR") +
  theme_minimal()

MSWEP = ggplot(region_mswep, aes(x = x, y = y, fill = precipitation)) +
  geom_raster() +
  scale_fill_viridis_c(option = "C",
                       limits = c(0, 1000),
                       breaks = c(0, 200, 400, 600, 800, 1000)) +
  coord_equal() +
  labs(subtitle = "MSWEP") +
  theme_minimal()

JAXA = ggplot(region_jaxa, aes(x = lon, y = lat, fill = precip)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C", 
                       name = "Precip (mm/month)",
                       limits = c(0, 1000),
                       breaks = c(0, 200, 400, 600, 800, 1000)) +
  coord_equal() +
  labs(subtitle = "JAXA GSMAP") +
  theme_minimal()

CMORPH = ggplot(oct16_sum_grid, aes(x = lon, y = lat, fill = precip)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C", 
                       name = "Precip (mm/month)",
                       limits = c(0, 1000),
                       breaks = c(0, 200, 400, 600, 800, 1000)) +
  coord_equal() +
  labs(subtitle = "CMORPH") +
  theme_minimal()

grid.arrange(NASA, ERA5, CDR, MSWEP, JAXA, CMORPH, nrow = 3)

