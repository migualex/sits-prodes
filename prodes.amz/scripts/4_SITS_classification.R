# ============================================================
#  Classification of Vector Data Cube
# ============================================================

## I. Load Required Libraries
library(sits)
library(tibble)
library(ggplot2)

## II. Define the date and time for the start of processing
date_process <- format(Sys.Date(), "/%Y_%m_%d_")
time_process <- format(Sys.time(), "%Hh%Mm_", tz = "America/Sao_Paulo")
process_version <- paste0(date_process, time_process)

## III. Define the paths for files and folders needed in the processing
model_name    <- "2025_12_12_08h46m_RF_1y_012015_012014_013015_013014_with_df_mask.rds" #add the model name
seg_version   <- "snic-1ymlme-rectangular-compactness-05"# SITS recognizes "underline" as a separator of information. Use only for this purpose.
vector_path   <- "data/segments"
class_path    <- "data/class"
rds_path      <- paste0("data/rds/model/random_forest/", model_name)
mixture_path  <- "data/raw/mixture_model"

var <- stringr::str_extract(basename(rds_path), "(with|no)-df-mask")

# ============================================================
# 1. Define and Load Data Cubes
# ============================================================

# Step 1.1 -- Create a classification cube from a collection
cube <- sits_cube(
  source      = "BDC",
  collection  = "SENTINEL-2-16D",
  bands       = c('B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12', 'NDVI', 'NBR', 'EVI', 'CLOUD'),
  tiles       = "012014",
  start_date  = "2024-08-01",
  end_date    = "2025-07-31",
  progress    = TRUE
)

# Step 1.1.1 -- Concatenates all the names of the classification tiles into a single string separated by '-'
tiles_class <- paste(cube$tile, collapse = "-")

datas <- sits_timeline(cube)
qtd_anos <- paste0(floor(lubridate::interval(datas[1], datas[length(datas)]) / lubridate::years(1)), "y")

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
version <- paste("rf", qtd_anos, tiles_class, var, sep = "-")

# Step 2.3 -- Classify segments according to the probabilities and calculate the process duration
sits_classify_start <- Sys.time()
class_prob <- sits_classify(
  data        = local_segs_cube,
  ml_model    = rf_model,
  multicores  = 4,
  memsize     = 80,
  output_dir =  class_path,
  version     = version,
  n_sam_pol   = 16, #  Number of time series per segment to be classified (integer, min = 10, max = 50)
  verbose     = TRUE,
  progress    = TRUE
)
sits_classify_end <- Sys.time()
temp_process_sits_classify <- round(sits_classify_end-sits_classify_start,2)
temp_process_sits_classify

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
  multicores  = 4,
  memsize     = 80,
  progress    = TRUE
)

print("Classification has finished")