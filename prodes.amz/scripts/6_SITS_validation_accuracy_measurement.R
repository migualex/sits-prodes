# ============================================================
#  Validation and accuracy measurement
# ============================================================

## I. Load Required Libraries
library(tibble)
library(sits)
library(terra)

## II. Define the date and time for the start of processing
date_process <- format(Sys.Date(), "/%Y_%m_%d_")
time_process <- format(Sys.time(), "%Hh%Mm_", tz = "America/Sao_Paulo")
process_version <- paste0(date_process, time_process)

## III. Define the paths for files and folders needed in the processing
data_dir <- "data/class"
val_samples_dir <- "data/raw/validation_samples"
train_samples_dir <- "data/raw/training_samples"
mixture_path  <- "data/raw/mixture_model"

tiles_class <- paste(cube$tile, collapse = "-")
dates <- sits_timeline(cube)
qtd_anos <- paste0(floor(lubridate::interval(dates[1], dates[length(dates)]) / lubridate::years(1)), "y")
version <- paste("rf", qtd_anos, tiles_class, var, sep = "-")

# ============================================================
# 1. Load probabilities cube
# ============================================================

# Step 1.1 -- Define the classes of the probability cube
labels <- c("DESMAT_CORTE_RASO" = "Clear_Cut_Bare_Soil",
            "DESMAT_CORTE_RASO_DM" = "Clear_Cut_Bare_Soil_mask",
            "DESMAT_CORTE_RASO_FOGO" = "Clear_Cut_Burned_Area",
            "DESMAT_DEGRAD_FOGO" = "Deforested_by_Fire",
            "DESMAT_ARVORE_REMANESCE" = "Clear_Cut_Bare_Soil_remaining trees", 
            "FLORESTA" = "Forest",
            "FLO_DEGRAD_FOGO" = "Degraded_Forest_Fire",
            "FLO_DEGRAD" = "Degraded_Forest",
            "DESMAT_VEG" = "Clear_Cut_Vegetation",
            "DESMAT_VEG_DM" = "Clear_Cut_Vegetation_mask",
            "AGUA" = "Water",
            "NF" = "Non_Forest_Vegetation",
            "ROCHA" = "Rock",
            "WETLANDS" = "Wetland")


cube <- sits_cube(
  source      = "BDC",
  collection  = "SENTINEL-2-16D",
  bands       = c('B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12', 'NDVI', 'NBR', 'EVI', 'CLOUD'),
  tiles       = "012014",
  start_date  = "2022-08-01",
  end_date    = "2025-07-31",
  progress    = TRUE
)

mm_cube <- sits_cube(
  source = "BDC",
  tiles = '012014',
  collection = "SENTINEL-2-16D",
  bands = c("SOIL", "VEG", "WATER"),
  data_dir = mixture_path,
  progress = TRUE
)

cube_merge_lsmm_class <- sits_merge(mm_cube, cube)

# Step 1.2 -- Create a probability data cube from a file 
class_datacube <- sits_cube(
  source = "BDC",
  collection = "SENTINEL-2-16D",
  tiles = "012014",
  raster_cube = cube_merge_lsmm_class,
  vector_dir = data_dir,
  vector_band = "class",
  version = version,
  label
  parse_info  = c("X1", "X2", "tile", "start_date", "end_date", "band", "version"),
  progress = TRUE
)

# Step 1.2.1 --  Plot the classification map
plot(class_datacube)


# ============================================================
# 2. Accuracy assessment of classified images
# ============================================================

# Step 2.1 -- Sampling design
sampling_design <- sits_sampling_design(
  cube = class_datacube,
  expected_ua = c(
    "Clear_Cut_Bare_Soil" = 0.75,
    "Clear_Cut_Bare_Soil_mask" = 0.70,
    "Deforested_by_Fire" = 0.70,
    "Clear_Cut_Bare_Soil_remaining trees" = 0.70, 
    "Clear_Cut_Burned_Area" = 0.70, 
    "Clear_Cut_Vegetation_mask" = 0.70, 
    "Degraded_Forest_Fire" = 0.70,
    "Degraded_Forest" = 0.70,
    "Non_Forest_Vegetation" = 0.70,
    "Rock" = 0.70,
    "Forest" = 0.75,  
    "Clear_Cut_Vegetation" = 0.70,  
    "Water" = 0.70, 
    "Wetland" = 0.70
  ),
  alloc_options = c(120, 100),
  std_err = 0.01,
  rare_class_prop = 0.1
)


# Step 2.2 -- Show sampling design
sampling_design

# Step 2.3 -- Generate stratified random samples
samples_sf <- sits_stratified_sampling(
  cube = class_datacube,
  sampling_design = sampling_design,
  alloc = "alloc_120",
  multicores = 4 # adapt to your computer CPU core availability
  )

# Step 2.4 -- Save samples in a shapefile
sf::st_write(samples_sf, 
             file.path(val_samples_dir, "validation_samples.shp"), 
             append = FALSE
             )

# Step 2.5 -- Get ground truth points
ground_truth <- system.file(
  "class/samples_validation.csv", package = "sitsdata"
  )

# Step 2.6 -- Calculate accuracy according to Olofsson's method
area_acc <- sits_accuracy(probs_datacube_class, 
                          validation = ground_truth,
                          multicores = 4 # adapt to your computer CPU core availability
                          )

# Step 2.7 -- Print the area estimated accuracy
area_acc

# Step 2.8 -- Show confusion matrix
area_acc$error_matrix