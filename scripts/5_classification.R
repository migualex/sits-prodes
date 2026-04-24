# ============================================================
#  Classification of Vector Data Cube
# ============================================================

# Load required libraries
library(sits)
library(tibble)
library(ggplot2)
library(terra)
library(RColorBrewer)

# Define the parameters: These are user-defined variables
model_name    <- "rf-model_4t_012014-012015-013014-013015_1y_2024-08-01_2025-07-31_all-samples-new-pol-avg-false_2026-04-15_12h01m.rds"
seg_version   <- "lsmm-snic-spac10-comp03-pad0-rectangular"# SITS recognizes "underline" as a separator of information. Use only for this purpose.
start_date    <- "2024-08-01"
end_date      <- "2025-07-31"
tile          <- "012014" # one tile per classification run

# File and folder paths 
models <- c("rf"   = "random_forest",
            "xgb"  = "xgboost",
            "ltae" = "ltae",
            "tcnn" = "temp_cnn",
            "rnet" = "res_net",
            "lstm" = "ltsm")
model_type    <- stringr::str_split_i(model_name, "-", 1)
model_path    <- file.path("data/rds/model", models[model_type], model_name)
vector_path   <- "data/segments"
class_path    <- "data/class"
mixture_path  <- "data/raw/mixture_model"

# Identifier to distinguish this model run from previous runs
var <- stringr::str_split_i(model_name, "_", 7)

# ============================================================
# 1. Define and Load Data Cubes
# ============================================================

# Step 1.1 -- Create a classification cube from a collection
cube <- sits_cube(
  source      = "BDC",
  collection  = "SENTINEL-2-16D",
  bands       = c('B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12', 'NDVI', 'NBR', 'EVI', 'CLOUD'),
  tiles       = tile,
  start_date  = start_date,
  end_date    = end_date,
  progress    = TRUE)

# Step 1.2 -- Extract tiles and duration from the cube (in years)
no.years <- paste0(floor(lubridate::year(end_date) - lubridate::year(start_date)), "y")

# Step 1.3 -- Retrieve Mixture Model Cube from a predefined repository
mm_cube <- sits_cube(
  source      = "BDC",
  collection  = "SENTINEL-2-16D",
  bands       = c("SOIL", "VEG", "WATER"),
  tiles       = tile,
  data_dir    = mixture_path,
  start_date  = start_date,
  end_date    = end_date,
  progress    = TRUE)

# Step 1.4 -- Merge the Classification Cube with Mixture Model Cube
cube_merge_lsmm_class <- sits_merge(mm_cube, cube)

# Step 1.5 -- Create a local segmented cube based on previous segmentation results
local_segs_cube <- sits_cube(
  source      = "BDC",
  collection  = "SENTINEL-2-16D",
  raster_cube = cube_merge_lsmm_class,
  vector_dir  = vector_path,
  vector_band = "segments",
  version     = seg_version, 
  parse_info  = c("satellite", "sensor","tile", "start_date", 
                  "end_date", "band", "version", "X1"))

# Step 1.6 -- Create output directory per tile
tile_period_dir <- file.path(class_path, tile, "original_class")
dir.create(tile_period_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. Probability and Classification Mapping
# ============================================================

# Step 2.1 -- Retrieve the trained model
model <- readRDS(model_path)

# Step 2.2 -- Define the version name of probability file
version <- paste(model_type, no.years, var, sep = "-")

# Step 2.3 -- Classify segments according to the probabilities and calculate the process duration
sits_classify_start <- Sys.time()
class_prob <- sits_classify(
  data        = local_segs_cube,
  ml_model    = model,
  multicores  = 28,  # adapt to your computer CPU core availability
  memsize     = 180, # adapt to your computer memory availability
  output_dir  = tile_period_dir,
  version     = version,
  n_sam_pol   = 16, #  Number of time series per segment to be classified (integer, min = 10, max = 50)
  verbose     = TRUE,
  progress    = TRUE
)
sits_classify_end <- Sys.time()
sits_classify_time <- as.numeric(sits_classify_end - sits_classify_start, units = "secs")
sprintf("SITS classify process duration (HH:MM): %02d:%02d",
        as.integer(sits_classify_time / 3600),
        as.integer((sits_classify_time %% 3600) / 60))

# Step 2.4 -- Reconstruct vector cube with classification probabilities 
vector_cube <- sits_cube(
  source      = "BDC",
  collection  = "SENTINEL-2-16D",
  raster_cube = cube_merge_lsmm_class,
  vector_dir  = tile_period_dir,
  vector_band = "probs",
  version     = version, # do not use underline character
  parse_info  = c("X1", "X2", "tile", "start_date", "end_date", "band", "version")
)

# Step 2.5 -- Generate Final Classified Map of Segments
class_map <- sits_label_classification(
  cube        = class_prob,
  output_dir  = tile_period_dir,
  version     = version,
  multicores  = 28,  # adapt to your computer CPU core availability
  memsize     = 180, # adapt to your computer memory availability
  progress    = TRUE
)
print("Classification finished!")

# ============================================================
# 3. Uncertainty
# ============================================================

# Step 3.1 -- Define function to calculate entropy, rasterize and exclude .gpkg
compute_uncertainty_raster <- function(
  vector_cube,
  tile_period_dir,
  version,
  multicores = 28,
  memsize    = 180,
  delete_gpkg = TRUE
) {
  
  # Calculate uncertainty vector cube
  uncertainty <- sits_uncertainty(
    vector_cube,
    type       = "entropy",
    multicores = multicores,
    memsize    = memsize,
    output_dir = tile_period_dir,
    version    = version,
    progress   = TRUE
  )
  
  # List entropy .gpkg files and get the most recent one
  uncertainty_files <- list.files(
    path      = tile_period_dir,
    pattern   = "entropy.*\\.gpkg$",
    full.names = TRUE
  )
  uncertainty_file <- uncertainty_files[which.max(file.info(uncertainty_files)$mtime)]
  
  # Read the segment polygons file with entropy
  uncertainty_polygons <- terra::vect(uncertainty_file)
  
  # Create a raster template based on uncertainty_polygons
  raster_template <- terra::rast(
    terra::ext(uncertainty_polygons),
    res = terra::res(terra::rast(vector_cube$file_info[[1]]$path[1])),
    crs = terra::crs(uncertainty_polygons)
  )
  
  # Rasterize entropy values and scale to UINT16
  uncertainty_raster       <- terra::rasterize(uncertainty_polygons, raster_template, field = "entropy")
  uncertainty_raster_uint16 <- round(uncertainty_raster * 10000)
  
  # Plot
  plot(
    uncertainty_raster_uint16,
    col     = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))(100),
    maxcell = terra::ncell(uncertainty_raster_uint16),
    main    = "Uncertainty Map - Full Resolution"
  )
  
  # Save as .tif (UINT16, DEFLATE compressed)
  tile_period_dir <- file.path(class_path, tile, "original_class")
  
  tif_path <- file.path(
    tile_period_dir,
    paste0(tools::file_path_sans_ext(basename(uncertainty_file)), ".tif")
  )
  
  terra::writeRaster(
    uncertainty_raster_uint16,
    filename  = tif_path,
    datatype  = "INT2U",
    overwrite = TRUE,
    gdal      = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9"),
    progress  = TRUE
  )
  
  # Delete the source .gpkg files (optional)
  if (delete_gpkg) {
    removed <- file.remove(uncertainty_files)
    message("Deleted .gpkg files: ", paste(uncertainty_files[removed], collapse = ", "))
  }
  
  message("Uncertainty raster saved: ", tif_path)
  invisible(tif_path)
}

# Step 3.1 -- Run function to calculate entropy, rasterize and exclude .gpkg
compute_uncertainty_raster(
  vector_cube     = vector_cube,
  tile_period_dir = class_path,
  version         = version,
  multicores      = 28, # adapt to your computer CPU core availability
  memsize         = 180, # adapt to your computer CPU core availability
  delete_gpkg     = FALSE  # Keep the .gpkg file if you want to inspect it beforehand
)
