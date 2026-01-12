## I. Load Required Libraries
library(tibble)
library(sits)
library(sitsdata)

## II. Set tempdir 
tempdir_r <- "~/Gustavo/SITS/sits-prodes/prodes.amz/data"
dir.create(tempdir_r, showWarnings = FALSE)

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

# Step 1.2 -- Directory where the data is stored 
data_dir <- system.file("extdata/Rondonia-Class-2022-Mosaic/", package = "sitsdata")

# Step 1.3 -- Create a probability data cube from a file 
probs_datacube_class <- sits_cube(
  source = "MPC",
  collection = "SENTINEL-2-L2A",
  data_dir = data_dir,
  bands = "class",
  labels = labels,
  version = "mosaic"
)

# Step 1.3.1 --  Plot the classification map
plot(probs_datacube_class)

# ============================================================
# 2. Sampling design
# ============================================================

# Step 2.1 -- 
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

# ============================================================
# 3. Stratified random sampling
# ============================================================

# 3.1 -- Generate stratified samples
samples_sf <- sits_stratified_sampling(
  cube = probs_datacube_class,
  sampling_design = sampling_design,
  alloc = "alloc_120",
  multicores = 4
)

# 3.2 -- Save samples in a shapefile
sf::st_write(samples_sf, 
             file.path(tempdir_r, "samples.shp"), 
             append = FALSE
)

# ============================================================
# 4. Accuracy assessment of classified images
# ============================================================

# 4.1 -- Get ground truth points
valid_csv <- system.file(
  "class/samples_validation.csv", package = "sitsdata"
)

# 4.2 -- Calculate accuracy according to Olofsson's method
area_acc <- sits_accuracy(probs_datacube_class, 
                          validation = valid_csv,
                          multicores = 4)

# 4.2.1 -- Print the area estimated accuracy
area_acc

# 4.3 -- Show confusion matrix
area_acc$error_matrix