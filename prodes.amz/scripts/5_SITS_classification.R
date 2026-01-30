# ============================================================
#  Classification of Vector Data Cube
# ============================================================

## I. Load Required Libraries
library(sits)
library(tibble)
library(ggplot2)
library(terra)
library(RColorBrewer)

## II. Define the date and time for the start of processing
date_process <- format(Sys.Date(), "/%Y_%m_%d_")
time_process <- format(Sys.time(), "%Hh%Mm_", tz = "America/Sao_Paulo")
process_version <- paste0(date_process, time_process)

## III. Define the paths for files and folders needed in the processing
model_name    <- "RF-model_4-tiles-012015-012014-013015-013014_3y-period-2022-07-28_2025-07-28_with-df-mask-with-all-samples_2026-01-22_09h35m.rds" #add the model name
seg_version   <- "snic-1ymlme-rectangular-compactness-05"# SITS recognizes "underline" as a separator of information. Use only for this purpose.
vector_path   <- "data/segments"
class_path    <- "data/class"
rds_path      <- paste0("data/rds/model/random_forest/", model_name)
mixture_path  <- "data/raw/mixture_model"

var <- "with-df-mask"

# ============================================================
# 1. Define and Load Data Cubes
# ============================================================

# Step 1.1 -- Create a classification cube from a collection
cube <- sits_cube(
  source      = "BDC",
  collection  = "SENTINEL-2-16D",
  bands       = c('B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12', 'NDVI', 'NBR', 'EVI', 'CLOUD'),
  tiles       = "012014",
  start_date  = "2022-08-01",
  end_date    = "2025-07-31",
  progress    = TRUE
)

# Step 1.1.1 -- Concatenates all the names of the classification tiles into a single string separated by '-'
tiles_class <- paste(cube$tile, collapse = "-")

dates <- sits_timeline(cube)
no.years <- paste0(floor(lubridate::interval(dates[1], dates[length(dates)]) / lubridate::years(1)), "y")

# Step 1.2 -- Retrieve Mixture Model Cube from a predefined repository
mm_cube <- sits_cube(
  source = "BDC",
  tiles = '012014',
  collection = "SENTINEL-2-16D",
  bands = c("SOIL", "VEG", "WATER"),
  data_dir = mixture_path,
  progress = TRUE
)

# Step 1.3 -- Merge the Classification Cube with Mixture Model Cube
cube_merge_lsmm_class <- sits_merge(mm_cube, cube)

# Step 1.4 -- Create a local segmented cube based on previous segmentation results
local_segs_cube <- sits_cube(
  source      = "BDC",
  collection  = "SENTINEL-2-16D",
  raster_cube = cube_merge_lsmm_class,
  vector_dir  = vector_path,
  vector_band = "segments",
  version     = seg_version, 
  parse_info  = c("satellite", "sensor","tile", "start_date", "end_date", "band", "version")
)

# ============================================================
# 2. Probability and Classification Mapping
# ============================================================

# Step 2.1 -- Retrieve the trained model
rf_model <- readRDS(rds_path)

# Step 2.2 -- Define the version name of probability file
version <- paste("rf", no.years, tiles_class, var, sep = "-")

# Step 2.3 -- Classify segments according to the probabilities and calculate the process duration
sits_classify_start <- Sys.time()
class_prob <- sits_classify(
  data        = local_segs_cube,
  ml_model    = rf_model,
  multicores  = 10,  # adapt to your computer CPU core availability
  memsize     = 80, # adapt to your computer memory availability
  output_dir =  class_path,
  version     = version,
  n_sam_pol   = 16, #  Number of time series per segment to be classified (integer, min = 10, max = 50)
  verbose     = TRUE,
  progress    = TRUE
)
sits_classify_end <- Sys.time()
sits_classify_time <- as.numeric(sits_classify_end - sits_classify_start, units = "secs")
sprintf("SITS classify process duration (HH:MM): %02d:%02d", as.integer(sits_classify_time / 3600), as.integer((sits_classify_time %% 3600) / 60))

# Step 2.4 -- Reconstruct vector cube with classification probabilities 
vector_cube <- sits_cube(
  source      = "BDC",
  collection  = "SENTINEL-2-16D",
  raster_cube = cube_merge_lsmm_class,
  vector_dir  = class_path,
  vector_band = "probs",
  version     = version, # do not use underline character
  parse_info  = c("X1", "X2", "tile", "start_date", "end_date", "band", "version")
)

# Step 2.5 -- Generate Final Classified Map of Segments
class_map <- sits_label_classification(
  cube        = class_prob,
  output_dir  = class_path,
  version     = version,
  multicores  = 10,  # adapt to your computer CPU core availability
  memsize     = 80, # adapt to your computer memory availability
  progress    = TRUE
)

print("Classification has finished")

# ============================================================
# 3. Uncertainty
# ============================================================

# Step 3.1 -- Calculate uncertainty vector cube
uncertainty <- sits_uncertainty(
  vector_cube,
  type = "entropy",
  multicores = 10, # adapt to your computer CPU core availability
  memsize = 80, # adapt to your computer memory availability
  output_dir = class_path,
  version = version,
  progress = TRUE
)

# Step 3.2.1 -- List the paths of the '.gpkg' files in 'class_path' containing 'entropy'
uncertainty_files <- list.files(
  path = class_path, 
  pattern = "entropy.*\\.gpkg$", 
  full.names = TRUE
)

# Step 3.2.2 -- Sort by modification date and get the last one (most recent)
uncertainty_file <- uncertainty_files[which.max(file.info(uncertainty_files)$mtime)]

# Step 3.2.3 -- Read the segment polygons file with entropy
uncertainty_polygons <- vect(uncertainty_file)

# Step 3.3 -- Create a raster template based on uncertainty_polygons
raster_template <- rast(
  ext(uncertainty_polygons), 
  res = res(rast(vector_cube$file_info[[1]]$path[1])), # Get the actual resolution of the cube
  crs = crs(uncertainty_polygons)
)

# Step 3.4 -- Rasterize the values of the 'entropy' variable in the uncertainty vector file
uncertainty_raster <- rasterize(uncertainty_polygons, raster_template, field = "entropy")

# Step 3.5.1 -- Multiply by 10,000 to maintain accuracy
uncertainty_raster_uint16 <- round(uncertainty_raster * 10000)

# Step 3.5.2 -- Plot the resulting uncertainty cube
plot(uncertainty_raster_uint16, 
     col = colorRampPalette(rev(brewer.pal(11, "Spectral")))(100),
     maxcell = ncell(uncertainty_raster_uint16), # Usa TODOS os pixels
     main = "Uncertainty Map - Full Resolution")

# Step 3.6 -- Save the final file with the desired data type
writeRaster(
  uncertainty_raster_uint16, 
  filename = file.path(class_path, paste0(tools::file_path_sans_ext(basename(uncertainty_file)), "_raster.tif")),
  datatype = "INT2U",  # This is the code for Uint16
  overwrite = TRUE,
  gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9") # Additional compression to reduce file size
)