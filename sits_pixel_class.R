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
cube <- sits_cube(
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
cube_select <- sits_select(
  cube,
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
samples_rondonia_2025 <- sits_get_data(
  cube        = cube_select,
  samples     = samples_sf,
  start_date  = "2024-07-01",
  end_date    = "2025-08-20",
  label       = "label",
  multicores  = 50,
  progress    = TRUE
)

# Visualization of the temporal patterns of the classes
samples_rondonia_2025 |> 
  sits_select(bands = c("NDVI", "EVI"), start_date = '2024-07-01', end_date = '2025-08-20') |> 
  sits_patterns() |> 
  plot()

# ============================================================
# 5. Quality assessment using Self-Organizing Maps (SOM)
# ============================================================
# Cluster samples using SOM
som_cluster <- sits_som_map(
  samples_rondonia_2025,
  grid_xdim  = 12,
  grid_ydim  = 12,
  rlen       = 100,
  distance   = "dtw",
  som_radius = 2,
  mode       = "online"
)

# Visualize SOM clustering results
plot(som_cluster)

# Evaluate SOM cluster quality
som_eval <- sits_som_evaluate_cluster(som_cluster)
plot(som_eval)
print(som_eval)

# Clean noisy or mixed samples
all_samples <- sits_som_clean_samples(
  som_map            = som_cluster, 
  prior_threshold    = 0.6,
  posterior_threshold= 0.6,
  keep               = c("clean", "analyze", "remove")
)
plot(all_samples)

# Keep only clean and analyzable samples
new_samples_v2 <- sits_som_clean_samples(
  som_map            = som_cluster, 
  prior_threshold    = 0.6,
  posterior_threshold= 0.6,
  keep               = c("clean", "analyze")
)
summary(new_samples_v2)

# Re-run SOM on cleaned samples
som_cluster_new <- sits_som_map(
  new_samples_v2,
  grid_xdim  = 12,
  grid_ydim  = 12,
  rlen       = 100,
  distance   = "dtw",
  som_radius = 2,
  mode       = "online"
)

# Evaluate cleaned sample set
som_eval_new <- sits_som_evaluate_cluster(som_cluster_new)
plot(som_eval_new)
print(som_eval_new)

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
  data        = cube_s1_s2, 
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
    "desmat_solo"   = 15.70,
    "desmat_veg"    = 14.68,
    "desmat_flo"    = 19.30,
    "agua"          = 12.09,
    "degrad_flo"    = 15.37,
    "degrad_fogo"   = 14.26,
    "floresta"      = 30.64,
    "nf"            = 17.37,
    "wetlands"      = 16.06
  ),
  multicores   = 4,
  memsize      = 50,
  data_dir     = class_path,
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
  samples     = samples_rondonia_2025,
  folds       = 5,
  ml_method   = sits_rfor(),
  multicores  = 4
)
print(rfor_validate_mt)
