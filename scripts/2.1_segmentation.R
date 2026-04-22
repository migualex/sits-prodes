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

start_date <- '2025-07-01'
end_date   <- '2025-07-31'
tiles      <- '014002'

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

mm_cube_fraction <- sits_merge(mm_cube, cube)

# ============================================================
# 2. Image Segmentation
# ============================================================

grid_seeding <- "rectangular"
spacing      <- 10
compactness  <- 0.3
padding      <- 0

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
  memsize    = 180,
  multicores = 28,
  output_dir = segments_path,
  version    = version
)

# ============================================================
# 3. Helper Functions
# ============================================================

calculate_bbox_area <- function(cube) {
  width    <- as.numeric(cube$xmax[[1]] - cube$xmin[[1]])
  height   <- as.numeric(cube$ymax[[1]] - cube$ymin[[1]])
  area_km2 <- (width * height) / 1e6
  return(round(area_km2))
}

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

end_date_fixed   <- as.Date("2025-07-28")
start_date_month <- as.Date(start_date, "%Y-%m-%d")

# ============================================================
# CLEAN BEFORE LOOP START
# ============================================================

seg_files_init <- list.files(
  segments_path,
  pattern = paste0("_", tiles, "_.*\\.gpkg$"),
  full.names = TRUE,
  recursive = FALSE
)

if (length(seg_files_init) > 0) {
  file.remove(seg_files_init)
  cat("Initial segmentation files removed before loop.\n")
}

# ============================================================
# FIXED BBOX AREA (computed once)
# ============================================================

bbox_area <- calculate_bbox_area(cube)
cat(sprintf("BBOX area (fixed): %d km2\n", bbox_area))

prev_seg_file <- NULL

# ============================================================
# LOOP SAFETY LIMIT
# ============================================================

max_iter <- 5
iter <- 0

while (iter < max_iter) {
  
  iter <- iter + 1
  
  cat(sprintf("\nIteration %d / %d\n", iter, max_iter))
  
  cat(sprintf("Period    : %s to %s\n", start_date_month, end_date_fixed))
  cat(sprintf("BBOX area : %d km2\n", bbox_area))
  
  start_date_month <- floor_date(start_date_month - days(1), "month")
  
  cat(sprintf("Processing period: %s to %s\n",
              start_date_month, end_date_fixed))
  
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
  
  version <- paste0(
    "LSMM-SNIC-spac", spacing,
    "-comp", gsub("[.]", "", as.character(compactness)),
    "-pad", padding,
    "-", grid_seeding,
    "-", date_process
  )
  
  sits_segment(
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
  
  new_seg_file <- get_latest_segment_file(segments_path, tiles)
  
  if (is.null(new_seg_file)) {
    cat("Warning: no segment file found. Stopping.\n")
    break
  }
  
  cat(sprintf("Generated file: %s\n", basename(new_seg_file)))
  
  if (!is.null(prev_seg_file) &&
      file.exists(prev_seg_file) &&
      normalizePath(prev_seg_file) != normalizePath(new_seg_file)) {
    file.remove(prev_seg_file)
    cat(sprintf("Removed previous file: %s\n", basename(prev_seg_file)))
  }
  
  total_segment_area <- calculate_segments_area(new_seg_file)
  cat(sprintf("Segments area: %d km2\n", total_segment_area))
  
  if (bbox_area <= total_segment_area) {
    if (!is.null(prev_seg_file) && file.exists(prev_seg_file)) {
      file.remove(prev_seg_file)
      cat(sprintf("Removed incomplete file: %s\n",
                  basename(prev_seg_file)))
    }
    cat("Full coverage reached. Stopping loop.\n")
    break
  }
  
  prev_seg_file <- new_seg_file
}

if (iter >= max_iter && bbox_area > total_segment_area) {
  cat("Max iterations reached (5) without full coverage. Processing stopped.\n")
}
print("Segmentation complete.")