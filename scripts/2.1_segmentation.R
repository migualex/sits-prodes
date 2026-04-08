# ============================================================
# Segmentation Using Linear Spectral Mixture Model and Sentinel-2 Bands
# ============================================================

# Load Required Libraries
library(sf, lib.loc = "/opt/r/R/x86_64-pc-linux-gnu-library/4.4")
library(lubridate)
library(sits)

# Define the paths for files and folders needed in the processing
mixture_path  <- "data/raw/mixture_model"
segments_path <- "data/segments"

# Define the date and time for the start of processing
date_process  <- format(Sys.Date(), "%Y-%m-%d")

# ============================================================
# 1. Load Mixture Model Data Cubes from a collection
# ============================================================

# Step 1.1 -- Define start_date, end_date and tiles for cube
start_date <- '2025-07-01'
end_date   <- '2025-07-31'
tiles      <- '014002'

cube <- sits_cube(
  source     = "BDC",
  collection = "SENTINEL-2-16D",
  bands      = c('B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12', 'CLOUD'),
  tiles      = tiles,
  start_date = start_date,
  end_date   = end_date,
  progress   = TRUE)

# Step 1.2 -- Create a cube from fraction image features
mm_cube <- sits_cube(
  source     = "BDC",
  collection = "SENTINEL-2-16D",
  bands      = c("SOIL", "VEG", "WATER"),
  tiles      = tiles,
  data_dir   = mixture_path,
  start_date = start_date,
  end_date   = end_date,
  progress   = TRUE)

mm_cube_fraction <- sits_merge(mm_cube, cube)

# ============================================================
# 2. Segment the merged cube 
# ============================================================

# Define parameters for SNIC segmentation
grid_seeding <- "rectangular"
spacing      <- 10
compactness  <- 0.3
padding      <- 0

# Step 2.1 -- Define the version of your segmentation file
version <- paste0("LSMM-SNIC-spac", spacing, "-comp", gsub("[.]", "", as.character(compactness)), "-pad", padding, "-", grid_seeding, "-", date_process)

# Step 2.2 -- Segment sits cube using SNIC algorithm
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
  version    = version)

# ============================================================
# 3. Calculate the area of a cube in km2
# ============================================================

# Step 3.1 -- Calculate the bounding box area
calculate_bbox_area <- function(cube) {
  width    <- as.numeric(cube$xmax[[1]] - cube$xmin[[1]])
  height   <- as.numeric(cube$ymax[[1]] - cube$ymin[[1]])
  area_km2 <- (width * height) / 1e6
  return(round(area_km2))
}

# Step 3.2 -- Calculate the total area covered by segment files for a given tile
calculate_segments_area <- function(segments_path, tile_number) {
  # List segment files matching the tile number
  seg_files <- list.files(segments_path,
                          pattern    = paste0(tile_number, ".*\\.(gpkg|shp)$"),
                          full.names = TRUE,
                          recursive  = TRUE)
  
  # Return 0 if no files are found
  if (length(seg_files) == 0) return(0)
  
  # Read only files with valid geometries
  all_segments <- lapply(seg_files, function(f) {
    s <- sf::st_read(f, quiet = TRUE)
    if (nrow(s) > 0) return(s) else return(NULL)
  })
  
  all_segments <- do.call(rbind, all_segments)
  
  # Return 0 if no features are found
  if (is.null(all_segments) || nrow(all_segments) == 0) return(0)
  
  areas_m2  <- as.numeric(sf::st_area(all_segments))
  total_km2 <- sum(areas_m2) / 1e6
  return(round(total_km2))
}

# ============================================================
# 4. Loop to Expand the Period
# ============================================================

# Step 4.1 -- Define fixed end date and initialize start date for the loop
end_date_fixed   <- as.Date("2025-07-28")
start_date_month <- as.Date(start_date, "%Y-%m-%d")

while (TRUE) {
  
  # Step 4.2 -- Calculate and display BBOX area for the current period
  bbox_area <- calculate_bbox_area(cube)
  cat(sprintf("Period    : %s to %s\n", start_date_month, end_date_fixed))
  cat(sprintf("BBOX area : %d km2\n",   bbox_area))
  
  # Step 4.3 -- Move start date back by one month
  start_date_month <- floor_date(start_date_month - days(1), "month")
  cat(sprintf("Processing period: %s to %s\n", start_date_month, end_date_fixed))
  
  # Step 4.4 -- Update cubes for the new period
  new_cube <- sits_cube(
    source     = "BDC",
    collection = "SENTINEL-2-16D",
    bands      = c('B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12', 'CLOUD'),
    tiles      = tiles,
    start_date = format(start_date_month, "%Y-%m-%d"),
    end_date   = format(end_date_fixed, "%Y-%m-%d"),
    progress   = TRUE)
  
  mm_cube <- sits_cube(
    source     = "BDC",
    collection = "SENTINEL-2-16D",
    bands      = c("SOIL", "VEG", "WATER"),
    tiles      = tiles,
    data_dir   = mixture_path,
    start_date = format(start_date_month, "%Y-%m-%d"),
    end_date   = format(end_date_fixed, "%Y-%m-%d"),
    progress   = TRUE)
  
  mm_cube_fraction <- sits_merge(mm_cube, new_cube)
  
  # Step 4.5 -- Build version tag and segment the updated cube
  version <- paste0("LSMM-SNIC-spac", spacing, "-comp", gsub("[.]", "", as.character(compactness)), "-pad", padding, "-", grid_seeding, "-", format(start_date_month, "%Y%m"))
  
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
    version    = version)
  
  # Step 4.6 -- Remove outdated segmentation files for this tile
  all_files <- list.files(segments_path,
                          pattern    = paste0(tiles, ".*\\.(gpkg|shp)$"),
                          full.names = TRUE,
                          recursive  = TRUE)
  
  current_tag <- format(start_date_month, "%Y%m")
  old_files   <- all_files[!grepl(current_tag, all_files)]
  
  if (length(old_files) > 0) {
    file.remove(old_files)
  }
  
  # Step 4.7 -- Calculate total segmented area for this tile
  total_segment_area <- calculate_segments_area(segments_path, tiles)
  cat(sprintf("Segs area : %d km2\n", total_segment_area))
  
  # Step 4.8 -- Stop loop when full coverage is reached
  if (bbox_area <= total_segment_area) {
    cat("Full coverage reached. Stopping loop.\n")
    break
  }
}
print("Segmentation complete.")
