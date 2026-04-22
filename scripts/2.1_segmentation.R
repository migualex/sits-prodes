# ============================================================
# Segmentation Using Linear Spectral Mixture Model and Sentinel-2 Data
# ============================================================

# Load required libraries
library(sf, lib.loc = "/opt/r/R/x86_64-pc-linux-gnu-library/4.4")
library(lubridate)
library(sits)

# Define input/output paths
mixture_path  <- "data/raw/mixture_model"
segments_path <- "data/segments"

# Store processing date (used for versioning outputs)
date_process  <- format(Sys.Date(), "%Y-%m-%d")

# ============================================================
# 1. Load Data Cubes (Spectral Bands + Mixture Model Fractions)
# ============================================================

# Step 1.1 -- Define temporal range and tile ID
start_date <- '2025-07-01'
end_date   <- '2025-07-31'
tiles      <- '016003'

# Load Sentinel-2 spectral bands cube
cube <- sits_cube(
  source     = "BDC",
  collection = "SENTINEL-2-16D",
  bands      = c('B02', 'B03', 'B04', 'B05', 'B06', 'B07',
                 'B08', 'B8A', 'B11', 'B12', 'CLOUD'),
  tiles      = tiles,
  start_date = start_date,
  end_date   = end_date,
  progress   = TRUE
)

# Step 1.2 -- Load mixture model fraction cube (SOIL, VEG, WATER)
mm_cube <- sits_cube(
  source     = "BDC",
  collection = "SENTINEL-2-16D",
  bands      = c("SOIL", "VEG", "WATER"),
  tiles      = tiles,
  data_dir   = mixture_path,
  start_date = start_date,
  end_date   = end_date,
  progress   = TRUE
)

# Merge spectral bands with fraction data into a single cube
mm_cube_fraction <- sits_merge(mm_cube, cube)

# ============================================================
# 2. Image Segmentation (SNIC Algorithm)
# ============================================================

# Define SNIC segmentation parameters
grid_seeding <- "rectangular"
spacing      <- 10
compactness  <- 0.3
padding      <- 0

# Step 2.1 -- Build a version string for output files
version <- paste0(
  "LSMM-SNIC-spac", spacing,
  "-comp", gsub("[.]", "", as.character(compactness)),
  "-pad", padding,
  "-", grid_seeding,
  "-", date_process
)

# Step 2.2 -- Run segmentation on the merged cube
mm_cube_segments <- sits_segment(
  cube      = mm_cube_fraction,
  seg_fn    = sits_snic(
    grid_seeding = grid_seeding,
    spacing      = spacing,
    compactness  = compactness,
    padding      = padding
  ),
  memsize    = 180,
  multicores = 28,
  output_dir = segments_path,
  version    = version
)

# ============================================================
# 3. Helper Functions
# ============================================================

# Step 3.1 -- Compute bounding box area (km²)
calculate_bbox_area <- function(cube) {
  width    <- as.numeric(cube$xmax[[1]] - cube$xmin[[1]])
  height   <- as.numeric(cube$ymax[[1]] - cube$ymin[[1]])
  area_km2 <- (width * height) / 1e6
  return(round(area_km2))
}

# Step 3.2 -- Retrieve the most recent .gpkg file for a given tile
# Only matches exact tile ID (ignores subdirectories)
get_latest_segment_file <- function(segments_path, tile_number) {
  seg_files <- list.files(
    segments_path,
    pattern    = paste0("_", tile_number, "_.*\\.gpkg$"),
    full.names = TRUE,
    recursive  = FALSE
  )
  
  if (length(seg_files) == 0) return(NULL)
  
  seg_files[order(file.mtime(seg_files), decreasing = TRUE)][1]
}

# Step 3.3 -- Compute total segmented area (km²) from a .gpkg file
calculate_segments_area <- function(seg_file) {
  if (is.null(seg_file) || !file.exists(seg_file)) return(0)
  
  segments <- tryCatch({
    s <- sf::st_read(seg_file, quiet = TRUE)
    if (!is.null(s) && nrow(s) > 0) s else NULL
  }, error = function(e) {
    message(sprintf("Error reading file: %s\n  -> %s", seg_file, e$message))
    NULL
  })
  
  if (is.null(segments)) return(0)
  
  areas_m2  <- as.numeric(sf::st_area(segments))
  total_km2 <- sum(areas_m2, na.rm = TRUE) / 1e6
  return(round(total_km2))
}

# ============================================================
# 4. Iterative Loop to Expand Temporal Coverage
# ============================================================

# Step 4.1 -- Define fixed end date and initialize loop start date
end_date_fixed   <- as.Date("2025-07-28")
start_date_month <- as.Date(start_date, "%Y-%m-%d")

# Keeps track of the previous segmentation file
prev_seg_file <- NULL

while (TRUE) {
  
  # Step 4.2 -- Compute and display bounding box area
  bbox_area <- calculate_bbox_area(cube)
  cat(sprintf("Period    : %s to %s\n", start_date_month, end_date_fixed))
  cat(sprintf("BBOX area : %d km2\n", bbox_area))
  
  # Step 4.3 -- Move start date one month backward
  start_date_month <- floor_date(start_date_month - days(1), "month")
  cat(sprintf("Processing period: %s to %s\n",
              start_date_month, end_date_fixed))
  
  # Step 4.4 -- Reload cubes for the updated time range
  new_cube <- sits_cube(
    source     = "BDC",
    collection = "SENTINEL-2-16D",
    bands      = c('B02','B03','B04','B05','B06','B07',
                   'B08','B8A','B11','B12','CLOUD'),
    tiles      = tiles,
    start_date = format(start_date_month, "%Y-%m-%d"),
    end_date   = format(end_date_fixed, "%Y-%m-%d"),
    progress   = TRUE
  )
  
  mm_cube <- sits_cube(
    source     = "BDC",
    collection = "SENTINEL-2-16D",
    bands      = c("SOIL", "VEG", "WATER"),
    tiles      = tiles,
    data_dir   = mixture_path,
    start_date = format(start_date_month, "%Y-%m-%d"),
    end_date   = format(end_date_fixed, "%Y-%m-%d"),
    progress   = TRUE
  )
  
  mm_cube_fraction <- sits_merge(mm_cube, new_cube)
  
  # Step 4.5 -- Rebuild version tag and run segmentation
  version <- paste0(
    "LSMM-SNIC-spac", spacing,
    "-comp", gsub("[.]", "", as.character(compactness)),
    "-pad", padding,
    "-", grid_seeding,
    "-", date_process
  )
  
  mm_cube_segments <- sits_segment(
    cube      = mm_cube_fraction,
    seg_fn    = sits_snic(
      grid_seeding = grid_seeding,
      spacing      = spacing,
      compactness  = compactness,
      padding      = padding
    ),
    memsize    = 60,
    multicores = 14,
    output_dir = segments_path,
    version    = version
  )
  
  # Step 4.6 -- Get the latest segmentation file for the tile
  new_seg_file <- get_latest_segment_file(segments_path, tiles)
  
  if (is.null(new_seg_file)) {
    cat("Warning: no segment file found for tile", tiles,
        "after segmentation. Skipping iteration.\n")
    next
  }
  
  cat(sprintf("Generated file: %s\n", basename(new_seg_file)))
  
  # Step 4.7 -- Remove previous file once a new one is confirmed
  if (!is.null(prev_seg_file) &&
      file.exists(prev_seg_file) &&
      normalizePath(prev_seg_file) != normalizePath(new_seg_file)) {
    file.remove(prev_seg_file)
    cat(sprintf("Removed previous file: %s\n", basename(prev_seg_file)))
  }
  
  # Step 4.8 -- Compute segmented area
  total_segment_area <- calculate_segments_area(new_seg_file)
  cat(sprintf("Segments area: %d km2\n", total_segment_area))
  
  # Step 4.9 -- Stop when full coverage is achieved
  if (bbox_area <= total_segment_area) {
    if (!is.null(prev_seg_file) && file.exists(prev_seg_file)) {
      file.remove(prev_seg_file)
      cat(sprintf("Removed incomplete file: %s\n",
                  basename(prev_seg_file)))
    }
    cat("Full coverage reached. Stopping loop.\n")
    break
  }
  
  # Step 4.10 -- Store current file for next iteration cleanup
  prev_seg_file <- new_seg_file
}
print("Segmentation complete.")
