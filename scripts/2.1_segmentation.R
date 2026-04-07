# ============================================================
# Segmentation Using Linear Spectral Mixture Model and Sentinel-2 Bands
# ============================================================

# Load Required Libraries
library(sf)
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
tiles      <- '012014'

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
# 2. Segment the mixture model cube based on fraction image features
# ============================================================

# Step 2.1 -- Define parameters for SNIC segmentation
spacing      <- 10
compactness  <- 0.3
padding      <- 0
grid_seeding <- "rectangular"

# Step 2.2 -- Define the version of your segmentation file
version <- paste0("LSMM-SNIC-spac", spacing, "-comp", gsub("\\.", "", as.character(compactness)), "-pad", padding, "-", grid_seeding, "-", date_process)

# Step 2.3 -- Segment sits cube using SNIC algorithm
sits_segment_start <- Sys.time()
mm_cube_segments <- sits_segment(
  cube      = mm_cube_fraction,
  seg_fn    = sits_snic(
    grid_seeding = grid_seeding,
    spacing      = spacing,
    compactness  = compactness,
    padding      = padding
  ),
  memsize    = 180, # adapt to your computer memory availability
  multicores = 28,  # adapt to your computer CPU core availability
  output_dir = segments_path,
  version    = version
)
sits_segment_end  <- Sys.time()
sits_segment_time <- as.numeric(sits_segment_end - sits_segment_start, units = "secs")
sprintf("SITS segment process duration (HH:MM): %02d:%02d", as.integer(sits_segment_time / 3600), as.integer((sits_segment_time %% 3600) / 60))
print("Segmentation complete.")

# ============================================================
# 3. Calculate Area of BBOX and Compare with Segment Areas
# ============================================================

# Step 3.1 -- Define a function to calculate the area of BBOX
calculate_bbox_area <- function(tile) {
  bbox <- st_bbox(as_spatial(tile))
  width <- abs(bbox$xmax - bbox$xmin)
  height <- abs(bbox$ymax - bbox$ymin)
  return(width * height)
}

# Step 3.2 -- Initialize variables for the loop
start_date_month <- as.Date(start_date, "%Y-%m-%d")
end_date_month <- as.Date(end_date, "%Y-%m-%d")

while (TRUE) {
  # Calculate the BBOX area of the tile
  bbox_area <- calculate_bbox_area(cube)
  
  # Calculate the sum of areas of segments
  segment_areas <- sits_area(mm_cube_segments)
  total_segment_area <- sum(segment_areas)
  
  # Print the current month, BBOX area and total segment area for comparison
  cat(sprintf("Date Range: %s to %s\n", start_date_month, end_date_month))
  cat(sprintf("BBOX Area: %.2f\n", bbox_area))
  cat(sprintf("Total Segment Area: %.2f\n", total_segment_area))
  
  # Check if BBOX area is greater than total segment area
  if (bbox_area > total_segment_area) {
    # Increase the period to the next month
    start_date_month <- end_date_month + days(1)
    end_date_month <- end_of_month(start_date_month)
    
    # Update the cube and segmentation with the new date range
    cube <- sits_cube(
      source     = "BDC",
      collection = "SENTINEL-2-16D",
      bands      = c('B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12', 'CLOUD'),
      tiles      = tiles,
      start_date = format(start_date_month, "%Y-%m-%d"),
      end_date   = format(end_date_month, "%Y-%m-%d"),
      progress   = TRUE
    )
    
    mm_cube <- sits_cube(
      source     = "BDC",
      collection = "SENTINEL-2-16D",
      bands      = c("SOIL", "VEG", "WATER"),
      tiles      = tiles,
      data_dir   = mixture_path,
      start_date = format(start_date_month, "%Y-%m-%d"),
      end_date   = format(end_date_month, "%Y-%m-%d"),
      progress   = TRUE
    )
    
    mm_cube_fraction <- sits_merge(mm_cube, cube)
    
    mm_cube_segments <- sits_segment(
      cube      = mm_cube_fraction,
      seg_fn    = sits_snic(
        grid_seeding = grid_seeding,
        spacing      = spacing,
        compactness  = compactness,
        padding      = padding
      ),
      memsize    = 180, # adapt to your computer memory availability
      multicores = 28,  # adapt to your computer CPU core availability
      output_dir = segments_path,
      version    = version
    )
    
    cat("Updated cube and segmentation with new date range.\n")
  } else {
    break
  }
}

# Final print statement when the loop ends
cat("BBOX area is equal or greater than total segment area. Process finished.\n")
