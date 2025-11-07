# ============================================================
#  Object-based Time Series Image Analysis
# ============================================================
## 1. Load Required Libraries ----
library(sits)
library(sitsdata)
library(sf)
library(tibble)
library(dplyr)
library(rstac)

vector_path   <- "~/SITs/amazônia/segmentação"
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

# Create a local segmented cube based on previous segmentation results
local_segs_cube <- sits_cube(
  source      = "BDC",
  collection  = "SENTINEL-2-16D",
  raster_cube = cube_all,
  vector_dir  = vector_path,
  vector_band = "segments",
  version     = "-step10-comp-10-mlme",
  parse_info  = c("X1", "X2", "tile", "start_date", "end_date", "band", "version")
)

# Display time series timeline
sits_timeline(local_segs_cube)

# ============================================================
# 3. Load and Explore Sample Data
# ============================================================
# Load sample shapefile (reference samples)
samples_sf  <- st_read(sample_path)

# Basic info about the samples
cat("Total samples:", nrow(samples_sf), "\n")
cat("Dimensions (rows x columns):", dim(samples_sf), "\n\n")

# Distribution of samples by class
table(samples_sf$label)

# ============================================================
# 4. Extract Time Series from Samples
# ============================================================
samples_rondonia <- sits_get_data(
  cube        = local_segs_cube,
  samples     = samples_sf,
  start_date  = "2024-07-01",
  end_date    = "2025-08-20",
  label       = "label",
  n_sam_pol   = 16,       # number of polygons per sample
  multicores  = 50,       # parallel processing
  progress    = TRUE
)

# Visualization of the temporal patterns of the classes
samples_rondonia |> 
  sits_select(bands = c("NDVI", "EVI"), start_date = '2024-07-01', end_date = '2025-08-20') |> 
  sits_patterns() |> 
  plot()

# ============================================================
# 5. Train Classification Model
# ============================================================
set.seed(03022024)  # for reproducibility

# Train Random Forest model
rf_model <- sits_train(
  samples   = samples_rondonia,
  ml_method = sits_rfor()
)

# Plot variable importance
plot(rf_model)

# ============================================================
# 6. Classification and Probability Mapping
# ============================================================
# Apply model to cube (classification by object)
class_prob <- sits_classify(
  data        = cube_one, # only tile 12014
  ml_model    = rf_model,
  multicores  = 4,
  memsize     = 50,
  output_dir  = class_path,
  version     = "rev",
  progress    = TRUE
)

# Reconstruct vector cube with classification probabilities
vector_cube <- sits_cube(
  source      = "BDC",
  collection  = "SENTINEL-2-16D",
  raster_cube = cube_all, # all tiles
  vector_dir  = class_path,
  vector_band = "probs",
  parse_info  = c("X1", "X2", "tile", "start_date", "end_date", "band", "version"),
  version     = "rev"
)

# ============================================================
# 7. Generate Final Classified Map
# ============================================================
class_map <- sits_label_classification(
  cube        = class_prob,
  output_dir  = class_path,
  version     = "rev",
  multicores  = 4,
  memsize     = 50,
  progress    = TRUE
)