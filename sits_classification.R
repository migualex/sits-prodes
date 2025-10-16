# --- 1. Install libraries and packages ---

library(sits)
library(sitsdata)
library(sf)
library(tibble)
library(dplyr)
library(rstac)


# --- 2. Directories, folders and paths ---

images_path1 <- "amazônia/imagens/sentinel-1/"
images_path2 <- "amazônia/imagens/sentinel-2/"
samp_path <- "amazônia/amostras/Merge_dsm_Degrad_prog_Prodes2025_Mg_sitsbook.shp"
classification_path <- "amazônia/classificação/"


# --- 3. Cube ---

## Sentinel-2 BDC
roi <- c(lon_min = -64.2072683839457312, lat_min = -8.5591214574530330, 
         lon_max = -63.2039174048757957, lat_max = -7.5841897221631491)

cube <- sits_cube(
  source = "BDC",
  collection = "SENTINEL-2-16D",
  bands = c("B02", "B03", "B04", "B8A", "B11", "B12"),
  roi = roi,
  start_date = "2024-07-27",
  end_date = "2025-06-26",
  out_dir = images_path2
)

cube_select <- sits_select(
  cube,
  bands = c("B02", "B03", "B04", "B8A", "B11", "B12"),
  roi = roi,
  start_date = "2024-07-27",
  end_date = "2025-06-26",
  data_dir = images_path2
)

plot(cube_select, red = "B11", green = "B8A", blue = "B02", date = "2025-06-10")


# Sentinel-2 MPC
cube_s2 <- sits_cube(
  source = "MPC",
  collection = "SENTINEL-2-L2A",
  bands = c("B02", "B8A", "B11"),
  tiles = c("20MLS","20MMS","20LLR","20LMR"),
  start_date = "2024-07-27",
  end_date = "2025-06-26"
)

cube_s2_reg <- sits_regularize(
  cube = cube_s2,
  period = "P1M",
  res = 10,
  tiles = c("20MLS","20MMS","20LLR","20LMR"),
  memsize = 12,
  multicores = 6,
  output_dir = images_path2
)

sits_timeline(cube_s2_reg)
plot(cube_s2_reg, red = "B11", green = "B8A", blue = "B02", date = "2025-06-01")


# --- 3.1 Sentinel-1 ---

cube_s1_rtc <- sits_cube(
  source = "MPC",
  collection = "SENTINEL-1-RTC",  
  bands = c("VV", "VH"),
  orbit = "descending",
  tiles = c("20MLS","20MMS","20LLR","20LMR"),
  start_date = "2024-07-01",
  end_date = "2024-07-31"
)

cube_s1_reg <- sits_regularize(
  cube = cube_s1_rtc,
  period = "P1M",
  res = 10,
  tiles = c("20MLS","20MMS","20LLR","20LMR"),
  memsize = 12,
  multicores = 6,
  output_dir = images_path1
)

plot(cube_s1_reg, band = "VH", palette = "Greys", scale = 0.7)

cube_s1_s2 <- sits_merge(cube_s2_reg, cube_s1_reg)
plot(cube_s1_s2, red = "B11", green = "B8A", blue = "VH")


# --- 3.2 NDWI ---

cube_select2 <- sits_apply(
  cube_select,
  NDWI = (B03 - B11) / (B03 + B11),
  output_dir = images_path3,
  progress = TRUE
)

plot(cube_select2, band = "NDWI", palette = "Blues")

max_ndwi_cube <- sits_reduce(
  cube_select2,
  NDWIMAX = t_max(NDWI),
  output_dir = images_path3,
  multicores = 20,
  progress = TRUE
)

plot(max_ndwi_cube, band = "NDWIMAX")


# --- 3.3 Modelo Linear de Mistura Espectral ---

em <- tibble::tribble(
  ~class, ~B02, ~B03, ~B04, ~B8A, ~B11, ~B12,
  "forest", 200, 352, 189, 2800, 1340, 546,
  "soil", 400, 650, 700, 3600, 3500, 1800,
  "water", 700, 1100, 1400, 850, 40, 26
)

reg_cube <- sits_mixture_model(
  data = cube_select,
  endmembers = em,
  multicores = 4,
  memsize = 12,
  output_dir = "atual/s2_teste/"
)

plot(reg_cube, red = "SOIL", green = "B8A", blue = "B02", date = "2025-06-26")


# --- 3.4 Cube exploration ---

bands_cube_select <- sits_bands(cube_select)
cat("Selected bands:\n", bands_cube_select, "\n\n")

timeline <- sits_timeline(cube_select)
cat("Cube timeline:\n")
print(timeline)
cat("\n")

first_date <- timeline[1]
cat("First date of the cube:\n")
print(first_date)
cat("\n")

last_date <- timeline[length(timeline)]
cat("Last date of the cube:\n")
print(last_date)


# --- 4. Sample analysis ---

sample_path <- samp_path3
samples_sf <- st_read(sample_path)

n_samples <- nrow(samples_sf)
cat("Total number of samples:", n_samples, "\n\n")

dim_samples <- dim(samples_sf)
cat("Data dimensions (rows x columns):", dim_samples, "\n")

table(samples_sf$label)
sits_view(samples_sf)


# --- 5.1 Sample analysis - Extracting time series ---

samples_rondonia_2025 <- sits_get_data(
  cube = reg_cube,
  samples = samples_sf,
  start_date = "2024-07-27",
  end_date = "2025-06-26",
  label = "label",
  multicores = 20,
  progress = TRUE
)

saveRDS(samples_rondonia_2025, "atual/samples_rondonia.rds")

samples_teste <- readRDS("atual/samples_rondonia.rds")
plot(sits_patterns(samples_teste))


# --- 5.2 Visualization ---

plot(samples_rondonia_2025)


# --- 5.3 Sample quality assessment ---

som_cluster <- sits_som_map(
  samples_rondonia_2025,
  grid_xdim = 12,
  grid_ydim = 12,
  rlen = 100,
  distance = "dtw",
  som_radius = 2,
  mode = "online"
)

plot(som_cluster)

som_eval <- sits_som_evaluate_cluster(som_cluster)
som_eval
plot(som_eval)

all_samples <- sits_som_clean_samples(
  som_map = som_cluster, 
  prior_threshold = 0.6,
  posterior_threshold = 0.6,
  keep = c("clean", "analyze", "remove")
)

plot(all_samples)

new_samples_v2 <- sits_som_clean_samples(
  som_map = som_cluster, 
  prior_threshold = 0.6,
  posterior_threshold = 0.6,
  keep = c("clean", "analyze")
)

summary(new_samples_v2)

som_cluster_new <- sits_som_map(
  new_samples_v2,
  grid_xdim = 12,
  grid_ydim = 12,
  rlen = 100,
  distance = "dtw",
  som_radius = 2,
  mode = "online"
)

som_eval_new <- sits_som_evaluate_cluster(som_cluster_new)
som_eval_new
plot(som_eval_new)

balanced_samples <- sits_reduce_imbalance(
  samples = samples_rondonia_2025,
  n_samples_over = 100,
  n_samples_under = 150,
  multicores = 4
)

summary(balanced_samples)
balanced_samples


# --- 6. Classification models ---

set.seed(03022024)

rf_model <- sits_train(
  samples = balanced_samples, 
  ml_method = sits_rfor()
)

plot(rf_model)
saveRDS(rf_model, "atual/rf_model.rds")

rf_model_teste <- readRDS("atual/rf_model.rds")
plot(rf_model_teste)

class_prob <- sits_classify(
  data = cube_s1_s2, 
  ml_model = rf_model,
  multicores = 20,         
  memsize = 50,            
  output_dir = classification_path3,
  version = "vp",
  progress = TRUE
)

rondonia_smooth <- sits_smooth(
  cube = class_prob,
  multicores = 20,
  memsize = 50,
  version = 'vp',
  output_dir = classification_path3,
  progress = TRUE
)

sits_view(rondonia_smooth)

class_map <- sits_label_classification(
  cube = rondonia_smooth,
  memsize = 50,
  multicores = 20,
  version = "vp",
  output_dir = classification_path3,
  progress = TRUE
)

sits_view(class_map)


# --- 6.3 Uncertainty cube ---

s2_cube_uncert <- sits_uncertainty(
  cube = rondonia_smooth,
  type = "margin",
  memsize = 50,
  multicores = 20,
  version = "vp",
  output_dir = classification_path3,
  progress = TRUE
)

plot(s2_cube_uncert)


# --- K-fold validation ---

rfor_validate_mt <- sits_kfold_validate(
  samples = samples_rondonia_2025,
  folds = 5,
  ml_method = sits_rfor(),
  multicores = 5
)

rfor_validate_mt


# --- Active learning ---

new_samples <- sits_uncertainty_sampling(
  uncert_cube = s2_cube_uncert,
  n = 200,
  min_uncert = 0.5,
  sampling_window = 10
)

sits_view(new_samples)

new_samples_sf <- st_as_sf(
  new_samples,
  coords = c("longitude", "latitude"),
  crs = 4674
)

st_write(
  new_samples_sf,
  dsn = "atual/amostras/new_samples_uncertainty.shp",
  delete_layer = TRUE
)
