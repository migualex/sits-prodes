# ============================================================
#  Classification of Vector Data Cube
# ============================================================

# ============================================================
# 1. Libraries, paths and some initial parameters
# ============================================================

# Step 1.1 -- Load Required Libraries
library(sits)
library(tibble)
library(ggplot2)
library(terra)
library(RColorBrewer)

# Step 1.2 -- Define the date and time for the start of processing
date_process <- format(Sys.Date(), "/%Y_%m_%d_")
time_process <- format(Sys.time(), "%Hh%Mm_", tz = "America/Sao_Paulo")
process_version <- paste0(date_process, time_process)

# Step 1.3 -- Define the paths for files and folders needed in the processing
model_name    <- "RF-model_2-tiles-014002-015002_0y-period-2024-07-27_2025-07-28_nf_2026-02-20_15h02m" #add the model name
seg_version   <- "lsmm-snic-spac10-comp05-pad0-hexagonal-2026-02-11"# SITS recognizes "underline" as a separator of information. Use only for this purpose.
vector_path   <- "data/segments"
class_path    <- "data/class"
rds_path      <- paste0("data/rds/model/random_forest/", model_name)
mixture_path  <- "data/raw/mixture_model"

# Step 1.4 -- Identifier to distinguish the file from previous versions 
var <- "nf"

# Step 1.5 -- Define time range
start_date    <- "2024-08-01"
end_date      <- "2025-07-31"


# ============================================================
# 2. Define and Load Data Cubes
# ============================================================

# Step 2.1 -- Create a classification cube from a collection
cube <- sits_cube(
  source      = "BDC",
  collection  = "SENTINEL-2-16D",
  bands       = c('B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12', 'NDVI', 'NBR', 'EVI', 'CLOUD'),
  tiles       = "014002",
  start_date  = start_date,
  end_date    = end_date,
  progress    = TRUE
)

# Step 2.2 -- Extract tiles, timeline and duration from the cube (in years)
tiles_class <- paste(cube$tile, collapse = "-")
dates <- sits_timeline(cube)
no.years <- paste0(floor(lubridate::interval(start_date, end_date) / lubridate::years(1)), "y")

# Step 2.3 -- Retrieve Mixture Model Cube from a predefined repository
mm_cube <- sits_cube(
  source = "BDC",
  tiles = '014002',
  collection = "SENTINEL-2-16D",
  bands = c("SOIL", "VEG", "WATER"),
  data_dir = mixture_path,
  start_date  = start_date,
  end_date    = end_date,
  progress = TRUE
)

# Step 2.4 -- Merge the Classification Cube with Mixture Model Cube
cube_merge_lsmm_class <- sits_merge(mm_cube, cube)

# Step 2.5 -- Create a local segmented cube based on previous segmentation results
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
# 3. Probability and Classification Mapping
# ============================================================

# Step 3.1 -- Retrieve the trained model
rf_model <- readRDS(rds_path)

# Step 3.2 -- Define the version name of probability file
version <- paste("rf", no.years, tiles_class, var, sep = "-")

# Step 3.3 -- Classify segments according to the probabilities and calculate the process duration
sits_classify_start <- Sys.time()
class_prob <- sits_classify(
  data        = local_segs_cube,
  ml_model    = rf_model,
  multicores  = 28,  # adapt to your computer CPU core availability
  memsize     = 180, # adapt to your computer memory availability
  output_dir =  class_path,
  version     = version,
  n_sam_pol   = 16, #  Number of time series per segment to be classified (integer, min = 10, max = 50)
  verbose     = TRUE,
  progress    = TRUE
)
sits_classify_end <- Sys.time()
sits_classify_time <- as.numeric(sits_classify_end - sits_classify_start, units = "secs")
sprintf("SITS classify process duration (HH:MM): %02d:%02d", as.integer(sits_classify_time / 3600), as.integer((sits_classify_time %% 3600) / 60))

# Step 3.4 -- Reconstruct vector cube with classification probabilities 
vector_cube <- sits_cube(
  source      = "BDC",
  collection  = "SENTINEL-2-16D",
  raster_cube = cube_merge_lsmm_class,
  vector_dir  = class_path,
  vector_band = "probs",
  version     = version, # do not use underline character
  parse_info  = c("X1", "X2", "tile", "start_date", "end_date", "band", "version")
)

# Step 3.5 -- Generate Final Classified Map of Segments
class_map <- sits_label_classification(
  cube        = class_prob,
  output_dir  = class_path,
  version     = version,
  multicores  = 28,  # adapt to your computer CPU core availability
  memsize     = 180, # adapt to your computer memory availability
  progress    = TRUE
)

print("Classification has finished")

# ============================================================
# 4. Uncertainty
# ============================================================

# Step 4.1 -- Calculate uncertainty vector cube
uncertainty <- sits_uncertainty(
  vector_cube,
  type = "entropy",
  multicores = 28, # adapt to your computer CPU core availability
  memsize = 180, # adapt to your computer memory availability
  output_dir = class_path,
  version = version,
  progress = TRUE
)

# Step 4.2 -- List the paths of the '.gpkg' files in 'class_path' containing 'entropy'
uncertainty_files <- list.files(
  path = class_path, 
  pattern = "entropy.*\\.gpkg$", 
  full.names = TRUE
)

# Step 4.3 -- Sort by modification date and get the last one (most recent)
uncertainty_file <- uncertainty_files[which.max(file.info(uncertainty_files)$mtime)]

# Step 4.4 -- Read the segment polygons file with entropy
uncertainty_polygons <- vect(uncertainty_file)

# Step 4.5 -- Create a raster template based on uncertainty_polygons
raster_template <- rast(
  ext(uncertainty_polygons), 
  res = res(rast(vector_cube$file_info[[1]]$path[1])), # Get the actual resolution of the cube
  crs = crs(uncertainty_polygons)
)

# Step 4.6 -- Rasterize the values of the 'entropy' variable in the uncertainty vector file
uncertainty_raster <- rasterize(uncertainty_polygons, raster_template, field = "entropy")

# Step 4.6.1 -- Multiply by 10,000 to maintain accuracy
uncertainty_raster_uint16 <- round(uncertainty_raster * 10000)

# Step 4.6.2 -- Plot the resulting uncertainty cube
plot(uncertainty_raster_uint16, 
     col = colorRampPalette(rev(brewer.pal(11, "Spectral")))(100),
     maxcell = ncell(uncertainty_raster_uint16), # Usa TODOS os pixels
     main = "Uncertainty Map - Full Resolution")


# Step 4.7 -- Save the final file with the desired data type
writeRaster(
  uncertainty_raster_uint16, 
  filename = file.path(class_path, paste0(tools::file_path_sans_ext(basename(uncertainty_file)), "_raster.tif")),
  datatype = "INT2U",  # This is the code for Uint16
  overwrite = TRUE,
  gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9") # Additional compression to reduce file size
)