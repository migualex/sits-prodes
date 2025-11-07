# ============================================================
#  Classification based on Pixel
# ============================================================

## 1. Load Required Libraries ----
library(sits)       # Satellite Image Time Series package
library(sitsdata)   # Provides example datasets for testing SITS
library(sf)         # Spatial vector data manipulation (modern 'sp' replacement)
library(tibble)     # Modern, tidy tables
library(dplyr)      # Data manipulation (filter, mutate, summarize, etc.)
library(rstac)      # Access to STAC catalogs (BDC, AWS, etc.)

sample_path   <- "~/SITs/amazônia/amostras/amostras_refinadas_012014.shp"
class_path    <- "~/SITs/amazônia/classificação"

# ============================================================
# 2. Define and Load Data Cubes
# ============================================================
# Create a cube from the BDC (Brazil Data Cube) collection
cube_all <- sits_cube(
  source      = "BDC",
  collection  = "SENTINEL-2-16D",
  bands       = c("B02", "B03", "B04", "B08", "B11", "B12", 
                  "NDVI", "EVI", "NBR", "CLOUD"),
  tiles       = c("012014", "12015", "13014", "13015"),
  start_date  = "2024-07-01",
  end_date    = "2025-08-20",
  progress    = TRUE
)

# Select and filter data from cube
cube_select <- sits_select(
  cube_all,
  bands       = c("B02", "B03", "B04", "B08", "B11", "B12", 
                  "NDVI", "EVI", "NBR", "CLOUD"),
  tiles       = c("012014", "12015", "13014", "13015"),
  start_date  = "2024-07-01",
  end_date    = "2025-08-20",
  progress    = TRUE
)

# Create a cube from the BDC (Brazil Data Cube) collection
cube_one <- sits_cube(
  source      = "BDC",
  collection  = "SENTINEL-2-16D",
  bands       = c("B02", "B03", "B04", "B08", "B11", "B12", 
                  "NDVI", "EVI", "NBR", "CLOUD"),
  tiles       = "012014",
  start_date  = "2024-07-01",
  end_date    = "2025-08-20",
  progress    = TRUE
)

# Select and filter data from cube
cube_select2 <- sits_select(
  cube_one,
  bands       = c("B02", "B03", "B04", "B08", "B11", "B12", 
                  "NDVI", "EVI", "NBR", "CLOUD"),
  tiles       = "012014",
  start_date  = "2024-07-01",
  end_date    = "2025-08-20",
  progress    = TRUE
)

# Inspect available bands and visualize a composite
sits_bands(cube_select)
plot(cube_select, red = "B11", green = "B8A", blue = "B02", date = "2025-06-10")

# ============================================================
# 3. Load and inspect training samples
# ============================================================
samples_sf  <- st_read(sample_path) 

# Sample summary
cat("Total number of samples:", nrow(samples_sf), "\n\n")
cat("Data dimensions (rows x columns):", dim(samples_sf), "\n")

# Samples per class
table(samples_sf$label)

# Visualize sample distribution
sits_view(samples_sf)

# ============================================================
# 4. Sample analysis and preprocessing
# ============================================================
## Extract time series for the sample points
samples_rondonia <- sits_get_data(
  cube        = cube_select,
  samples     = samples_sf,
  start_date  = "2024-07-01",
  end_date    = "2025-08-20",
  label       = "label",
  multicores  = 50,
  progress    = TRUE
)

# Visualization of the temporal patterns of the classes
samples_rondonia |> 
  sits_select(bands = c("NDVI", "EVI"), start_date = '2024-07-11', end_date = '2025-08-13') |> 
  sits_patterns() |> 
  plot()

# ============================================================
# 6. Classification modeling
# ============================================================
set.seed(03022024)

# Train Random Forest model
rf_model <- sits_train(
  samples   = balanced_samples, 
  ml_method = sits_rfor()
)

# Visualize variable importance
plot(rf_model)

# ============================================================
# 7. Classification and post-processing
# ============================================================
# Generate probability cube
class_prob <- sits_classify(
  data        = cube_select2, # only tile 12014
  ml_model    = rf_model,
  multicores  = 4,         
  memsize     = 50,            
  output_dir  = class_path,
  version     = "rev",
  progress    = TRUE
)

# ============================================================
# 8. Variance and smoothing filters
# ============================================================
# Compute class variance
rondonia_var <- sits_variance(
  cube          = class_prob,
  window_size   = 7,
  neigh_fraction= 0.50,
  output_dir    = class_path,
  multicores    = 4,
  memsize       = 50,
  version       = "rev", 
  progress      = TRUE
)
summary(rondonia_var)

# Apply smoothing to reduce classification noise
smooth_rondonia <- sits_smooth(
  class_prob,
  window_size   = 7,
  neigh_fraction= 0.5,
  smoothness = c(
    "desmat_solo"   = 15.70, # 100% variance
    "desmat_veg"    = 14.68, # 100% variance
    "desmat_flo"    = 19.30, # 100% variance
    "agua"          = 12.09, # 100% variance
    "degrad_flo"    = 15.37, # 100% variance
    "degrad_fogo"   = 14.26, # 100% variance
    "floresta"      = 30.64, # 100% variance
    "nf"            = 17.37, # 100% variance
    "wetlands"      = 16.06  # 100% variance
  ),
  multicores   = 4,
  memsize      = 50,
  output_dir   = class_path,
  version      = "rev",
  progress     = TRUE
)

# ============================================================
# 9. Generate final classification map
# ============================================================
class_map <- sits_label_classification(
  cube        = smooth_rondonia,
  memsize     = 50,
  multicores  = 4,
  version     = "rev",
  output_dir  = class_path,
  progress    = TRUE
)
sits_view(class_map)

# ============================================================
# 10. Model validation
# ============================================================
rfor_validate_mt <- sits_kfold_validate(
  samples     = samples_rondonia,
  folds       = 5,
  ml_method   = sits_rfor(),
  multicores  = 4
)
print(rfor_validate_mt)
