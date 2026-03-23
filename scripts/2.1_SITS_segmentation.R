# ============================================================
#  Segmentation from Linear Spectral Mixture Model (LSMM) for features training and classification
# ============================================================

# Load Required Libraries
library(sits)

# Define the paths for files and folders needed in the processing
mixture_path  <- "data/raw/mixture_model"
segments_path <- "data/segments"

# Define the date and time for the start of processing
date_process  <- format(Sys.Date(), "%Y-%m-%d")

# ============================================================
# 1. Load Mixture Model Data Cubes from a collection
# ============================================================

# Step 1.1 -- Define Start_date, End_date and tiles for cube 
start_date <- '2025-07-01'
end_date   <- '2025-07-31'
tiles      <- '012014'

cube <- sits_cube(
  source = "BDC",
  collection = "SENTINEL-2-16D",
  bands = c('B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12', 'CLOUD'),
  tiles = tiles,
  start_date = start_date,
  end_date = end_date,
  progress = TRUE)

# Step 1.2 -- Create a cube from Fraction Images features 
mm_cube <- sits_cube(
  source = "BDC",
  collection = "SENTINEL-2-16D",
  bands = c("SOIL", "VEG", "WATER"),
  tiles = tiles,
  data_dir = mixture_path,
  start_date = start_date,
  end_date = end_date,
  progress = TRUE)

mm_cube_fraction <- sits_merge(mm_cube, cube)

# ============================================================
# 2. Segment the mixture model cube based on Fraction Images features
# ============================================================

# Step 2.1 -- Define parameters for SNIC segmentation
spacing      <- 10
compactness  <- 0.5
padding      <- 0
grid_seeding <- "rectangular"

# Step 2.2 -- Define the version of your segmentation file
version <- paste0("LSMM-SNIC-spac",spacing,"-comp",gsub("\\.", "", as.character(compactness)),"-pad",padding,"-",grid_seeding,"-",date_process)

# Step 2.3 -- Segment sits cube using SNIC algorithm
sits_segment_start <- Sys.time()
mm_cube_segments <- sits_segment(
  cube = mm_cube_fraction,
  seg_fn = sits_snic(
    grid_seeding = grid_seeding,
    spacing = spacing,
    compactness = compactness,
    padding = padding
  ),
  memsize = 180, #adapt to your computer memory availability
  multicores = 28, #adapt to your computer CPU core availability
  output_dir = segments_path,
  version = version
)
sits_segment_end <- Sys.time()
sits_segment_time <- as.numeric(sits_segment_end - sits_segment_start, units = "secs")
sprintf("SITS segment process duration (HH:MM): %02d:%02d", as.integer(sits_segment_time / 3600), as.integer((sits_segment_time %% 3600) / 60))

print("Segmentation has finished")