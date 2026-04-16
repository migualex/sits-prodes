# ============================================================
# Segmentation Using Linear Spectral Mixture Model and Sentinel-2 Bands
# ============================================================

library(sf, lib.loc = "/opt/r/R/x86_64-pc-linux-gnu-library/4.4")
library(lubridate)
library(sits)

mixture_path  <- "data/raw/mixture_model"
segments_path <- "data/segments"
date_process  <- format(Sys.Date(), "%Y-%m-%d")

# ============================================================
# 1. Load Mixture Model Data Cubes from a collection
# ============================================================

start_date <- '2025-07-01'
end_date   <- '2025-07-31'
tiles      <- '014001'

cube <- sits_cube(
  source     = "BDC",
  collection = "SENTINEL-2-16D",
  bands      = c('B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12', 'CLOUD'),
  tiles      = tiles,
  start_date = start_date,
  end_date   = end_date,
  progress   = TRUE)

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

grid_seeding <- "rectangular"
spacing      <- 10
compactness  <- 0.3
padding      <- 0

version <- paste0("LSMM-SNIC-spac", spacing, "-comp",
                  gsub("[.]", "", as.character(compactness)),
                  "-pad", padding, "-", grid_seeding, "-", date_process)

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

calculate_bbox_area <- function(cube) {
  width    <- as.numeric(cube$xmax[[1]] - cube$xmin[[1]])
  height   <- as.numeric(cube$ymax[[1]] - cube$ymin[[1]])
  area_km2 <- (width * height) / 1e6
  return(round(area_km2))
}

calculate_segments_area <- function(segments_path, tile_number) {
  seg_files <- list.files(segments_path,
                          pattern    = paste0(tile_number, ".*\\.(gpkg|shp)$"),
                          full.names = TRUE,
                          recursive  = TRUE)
  
  if (length(seg_files) == 0) return(0)
  
  all_segments <- lapply(seg_files, function(f) {
    s <- sf::st_read(f, quiet = TRUE)
    if (nrow(s) > 0) return(s) else return(NULL)
  })
  
  all_segments <- do.call(rbind, all_segments)
  
  if (is.null(all_segments) || nrow(all_segments) == 0) return(0)
  
  areas_m2  <- as.numeric(sf::st_area(all_segments))
  total_km2 <- sum(areas_m2) / 1e6
  return(round(total_km2))
}

# ============================================================
# 4. Loop to Expand the Period
# ============================================================

end_date_fixed   <- as.Date("2025-07-28")
start_date_month <- as.Date(start_date, "%Y-%m-%d")

# BBOX is constant for the tile -- calculated once before the loop
bbox_area <- calculate_bbox_area(cube)
cat(sprintf("BBOX area (reference): %d km2\n", bbox_area))

while (TRUE) {
  
  # Step 4.1 -- Check current coverage
  total_segment_area <- calculate_segments_area(segments_path, tiles)
  cat(sprintf("Period    : %s to %s\n", start_date_month, end_date_fixed))
  cat(sprintf("BBOX area : %d km2\n",   bbox_area))
  cat(sprintf("Segs area : %d km2\n",   total_segment_area))
  
  # Step 4.2 -- Stop criterion
  if (bbox_area <= total_segment_area) {
    cat("Full coverage reached. Stopping loop.\n")
    break
  }
  
  # Step 4.3 -- Move back one month
  start_date_month <- floor_date(start_date_month - days(1), "month")
  cat(sprintf("Processing period: %s to %s\n", start_date_month, end_date_fixed))
  
  # Step 4.4 -- Update cubes for the expanded period
  new_cube <- sits_cube(
    source     = "BDC",
    collection = "SENTINEL-2-16D",
    bands      = c('B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12', 'CLOUD'),
    tiles      = tiles,
    start_date = format(start_date_month, "%Y-%m-%d"),
    end_date   = format(end_date_fixed,   "%Y-%m-%d"),
    progress   = TRUE)
  
  mm_cube <- sits_cube(
    source     = "BDC",
    collection = "SENTINEL-2-16D",
    bands      = c("SOIL", "VEG", "WATER"),
    tiles      = tiles,
    data_dir   = mixture_path,
    start_date = format(start_date_month, "%Y-%m-%d"),
    end_date   = format(end_date_fixed,   "%Y-%m-%d"),
    progress   = TRUE)
  
  mm_cube_fraction <- sits_merge(mm_cube, new_cube)
  
  # Step 4.5 -- New version tag with the current period's start month
  version <- paste0("LSMM-SNIC-spac", spacing, "-comp",
                    gsub("[.]", "", as.character(compactness)),
                    "-pad", padding, "-", grid_seeding,
                    "-", format(start_date_month, "%Y%m"),
                    "-", date_process)
  
  # Step 4.6 -- Segment the expanded period
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
  
  # Step 4.7 -- Remove files from previous versions (shorter period, now superseded)
  all_files   <- list.files(segments_path,
                            pattern    = paste0(tiles, ".*\\.(gpkg|shp)$"),
                            full.names = TRUE,
                            recursive  = TRUE)
  
  current_tag <- format(start_date_month, "%Y%m")
  old_files   <- all_files[!grepl(current_tag, all_files)]
  
  if (length(old_files) > 0) {
    file.remove(old_files)
  }
}
print("Segmentation complete.")
