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

tiles_class <- paste(cube$tile, collapse = "-")
dates <- sits_timeline(cube)
qtd_anos <- paste0(floor(lubridate::interval(dates[1], dates[length(dates)]) / lubridate::years(1)), "y")
version <- paste("rf", qtd_anos, tiles_class, var, sep = "-")

# ============================================================
# 1. Load probabilities cube
# ============================================================

# Step 1.1 -- Define the classes of the probability cube
labels <- c("1" = "Clear_Cut_Bare_Soil",
            "2" = "Clear_Cut_Burned_Area", 
            "3" = "Mountainside_Forest", 
            "4" = "Forest",
            "5" = "Riparian_Forest",
            "6" = "Clear_Cut_Vegetation",
            "7" = "Water",
            "8" = "Seasonally_Flooded",
            "9" = "Wetland")

# Step 1.2 -- Create a probability data cube from a file 
probs_datacube_class <- sits_cube(
  source = "BDC",
  collection = "SENTINEL-2-16D",
  data_dir = data_dir,
  bands = "class",
  labels = labels,
  version = version
)

# Step 1.2.1 --  Plot the classification map
plot(probs_datacube_class)


# ============================================================
# 2. Accuracy assessment of classified images
# ============================================================

# Step 2.1 -- Sampling design
sampling_design <- sits_sampling_design(
  cube = probs_datacube_class,
  expected_ua = c(
    "Clear_Cut_Bare_Soil" = 0.75,
    "Clear_Cut_Burned_Area" = 0.70, 
    "Mountainside_Forest" = 0.70, 
    "Forest" = 0.75,  
    "Riparian_Forest" = 0.70, 
    "Clear_Cut_Vegetation" = 0.70,  
    "Water" = 0.70, 
    "Seasonally_Flooded" = 0.70, 
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
  cube = probs_datacube_class,
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