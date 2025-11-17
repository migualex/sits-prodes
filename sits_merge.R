# ============================================================
#  Merge Sentinel-1 (SAR) and Sentinel-2 Data Cubes
# ============================================================

## 1. Load Required Libraries ----
library(sits)
library(sitsdata)
library(sf)
library(tibble)
library(dplyr)
library(rstac)

cube_path    <- "~/SITs/amazônia/imagens/s1/"

# ============================================================
# 2. Define ROI (Region of Interest)
# ============================================================
# ROI: Tile 012014 (Rondônia)
roi <- c(
  lon_min = -64.20726838394573,
  lat_min = -8.55912145745303,
  lon_max = -63.20391740487580,
  lat_max = -7.58418972216315
)

# ============================================================
# 3. Load Sentinel-1 (SAR) Cubes
# ============================================================
# Two temporal segments are used to ensure continuous coverage
cube_s1_part1 <- sits_cube(
  source      = "MPC",
  collection  = "SENTINEL-1-RTC",
  bands       = c("VV", "VH"),
  start_date  = "2024-08-01",
  end_date    = "2025-03-01",
  orbit       = "descending",
  roi         = roi,
  multicores  = 4
)

cube_s1_part2 <- sits_cube(
  source      = "MPC",
  collection  = "SENTINEL-1-RTC",
  bands       = c("VV", "VH"),
  start_date  = "2025-03-15",
  end_date    = "2025-07-31",
  orbit       = "descending",
  roi         = roi,
  multicores  = 4
)

# Merge both SAR cubes to create a single continuous cube
cube_s1 <- sits_merge(cube_s1_part1, cube_s1_part2)

# ============================================================
# 4. Load Sentinel-2 Cube
# ============================================================
cube_s2 <- sits_cube(
  source      = "BDC",
  collection  = "SENTINEL-2-16D",
  bands       = c('B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12', 'NDVI', 'NBR', 'EVI', 'CLOUD'),
  tiles       = "012014",
  start_date  = "2024-07-01",
  end_date    = "2025-08-12",
  progress    = TRUE
)

# ============================================================
# 5. Regularize Sentinel-1 Cube to Match Sentinel-2 Temporal Grid
# ============================================================
cube_s1_reg <- sits_regularize(
  cube        = cube_s1,
  period      = "P16D",                # same temporal frequency as S2
  res         = 10,                    # spatial resolution in meters
  tiles       = cube_s2$tile,
  timeline    = sits_timeline(cube_s2),
  memsize     = 50,
  multicores  = 4,
  grid_system = "BDC_SM_V2",
  output_dir  = cube_path
)

# Quick visual check of SAR data
# plot(cube_s1_reg, band = "VH", palette = "Greys", scale = 0.7)

# ============================================================
# 6. Merge Sentinel-1 and Sentinel-2 Cubes
# ============================================================
cube_s1_s2 <- sits_merge(cube_s2, cube_s1_reg)

# Confirm timeline alignment
# sits_timeline(cube_s1_s2)

# ============================================================
# 7. Visualization
# ============================================================
# Plot RGB-like composite with both optical (B11, B8A) and SAR (VH) information
plot(cube_s1_s2, red = "B11", green = "B8A", blue = "VH")
