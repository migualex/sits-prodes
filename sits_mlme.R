# ============================================================
#  Spectral Mixture Modeling and Image Segmentation
# ============================================================

## 1. Load Required Libraries ----
library(sits)
library(sitsdata)
library(sf)
library(tibble)
library(dplyr)
library(rstac)

seg_path            <- "~/SITs/amazônia/imagens/segmentacao/"
cube_mixture_path   <- "~/SITs/amazônia/imagens/mlme/"

# ============================================================
# 2. Load Sentinel-2 Data Cube
# ============================================================
cube <- sits_cube(
  source      = "BDC",
  collection  = "SENTINEL-2-16D",
  bands       = c("B02", "B03", "B04", "B08", "B8A",
                  "B11", "B12", "NDVI", "EVI", "NBR"),
  tiles       = "012014",
  start_date  = "2024-07-01",
  end_date    = "2025-08-12",
  progress    = TRUE
)

# Select relevant portion of the cube (explicit for clarity)
cube_select <- sits_select(
  cube,
  bands       = c("B02", "B03", "B04", "B08", "B8A",
                  "B11", "B12", "NDVI", "EVI", "NBR"),
  tiles       = "012014",
  start_date  = "2024-07-01",
  end_date    = "2025-08-12",
  progress    = TRUE
)

# Optional: inspect temporal coverage
# sits_timeline(cube_select)

# ============================================================
# 3. Define Endmembers for Linear Spectral Mixture Model
# ============================================================
# Typical endmembers: soil, forest, and water (reflectance values)
endmembers <- tibble::tribble(
  ~class,   ~B03, ~B04, ~B08, ~B11,
  "soil",     1275, 1501, 3958, 4035,
  "forest",   620,  207,  4450, 3123,
  "water",    329,  72,    86,   72
)

# ============================================================
# 4. Apply Linear Mixture Model
# ============================================================
cube_mixture <- sits_mixture_model(
  data        = cube_select,
  endmembers  = endmembers,
  multicores  = 4,
  memsize     = 50,
  output_dir  = cube_mixture_path,
  progress    = TRUE
)

# ============================================================
# 5. Select Fractional Endmember Bands
# ============================================================
cube_mixture_select <- sits_select(
  cube_mixture,
  bands       = c("FOREST", "SOIL", "WATER"),
  tiles       = "012014",
  start_date  = "2024-07-01",
  end_date    = "2025-08-12",
  progress    = TRUE
)

# ============================================================
# 6. Visualization of Endmember Fractions
# ============================================================
# Plot FOREST fraction
plot(cube_mixture, band = "FOREST", palette = "Greens")

# Plot SOIL fraction
plot(cube_mixture, band = "SOIL", palette = "OrRd")

# ============================================================
# 7. Image Segmentation (SLIC)
# ============================================================
segments_12014 <- sits_segment(
  cube        = cube_mixture_select,
  output_dir  = seg_path,
  seg_fn      = sits_slic(
    step        = 10,        # superpixel step size
    compactness = 0.5,       # spatial vs spectral weight
    dist_fun    = "euclidean",
    avg_fun     = "median",
    iter        = 10,
    minarea     = 10
  )
)

# ============================================================
# 8. Visualize Segmented Result
# ============================================================
# Note: updated object name (corrected typo from 'segments_20LMR')
plot(
  segments_12014,
  red   = "B11",
  green = "B8A",
  blue  = "B02",
  date  = "2025-06-11"
)