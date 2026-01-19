# ============================================================
#  Samples quality assessment, filtering and balancing
# ============================================================

## I. Load Required Libraries
library(sits)
library(sf)
library(tibble)
library(dplyr)
library(ggplot2)

## II. Define the date and time for the start of processing
date_process <- format(Sys.Date(), "/%Y_%m_%d_")
time_process <- format(Sys.time(), "%Hh%Mm_", tz = "America/Sao_Paulo")
process_version <- paste0(date_process, time_process)

## III. Define the paths for files and folders needed in the processing
sample_path   <- "data/raw/samples/" #add the sample file to the path
rds_path      <- "data/rds"
mixture_path  <- "data/raw/mixture_model"
plots_path    <- "data/plots"

# IV. Identifier to distinguish this model run from previous versions
var <- "_(with|no)_df_mask"

## V. Define a list with preference colours for each class
my_colours <- c(
  "OOB"                       = "black",
  "AGUA"                      = "#191ad7",
  "DESMAT_ARVORE_REMANESCE"   = "#e56c35",
  "DESMAT_CORTE_RASO"         = "#f01304",
  "DESMAT_CORTE_RASO_DM"      = "#f39c12",
  "DESMAT_DEGRAD_FOGO"        = "#a42900",
  "DESMAT_VEG"                = "#24fc15",
  "DESMAT_VEG_DM"             = "#e6b0aa",
  "FLO_DEGRAD"                = "#fbf909",
  "FLO_DEGRAD_FOGO"           = "#d1b007",
  "FLORESTA"                  = "#1e2f09",
  "NF"                        = "#fb0e9f",
  "ROCHA"                     = "#562917",
  "WETLANDS"                  = "#b779c6" 
)

# ============================================================
# 1. Define and Load Data Cubes
# ============================================================

# Step 1.1 -- Create a training cube from a collection
cube <- sits_cube(
  source      = "BDC",
  collection  = "SENTINEL-2-16D",
  bands       = c('B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12', 'NDVI', 'NBR', 'EVI', 'CLOUD'),
  tiles       = c("012014","012015","013014","013015"),
  start_date  = "2023-08-01",
  end_date    = "2025-07-31",
  progress    = TRUE
)

# Step 1.2 -- Calculate the number of years in the training cube
cube_dates <- sits_timeline(cube)
no.years <- paste0(floor(lubridate::interval(cube_dates[1], cube_dates[length(cube_dates)]) / lubridate::years(1)), "y")

# Step 1.3 -- Concatenates all the names of the training tiles into a single string separated by '-'
tiles_train <- paste(cube$tile, collapse = "-")

# Step 1.4 -- Retrieve Mixture Model Cube from a predefined repository
mm_cube <- sits_cube(
  source = "BDC",
  tiles = c('012014', '012015', '013014', '013015'),
  collection = "SENTINEL-2-16D",
  bands = c("SOIL", "VEG", "WATER"),
  data_dir = mixture_path,
  progress = TRUE
)

# Step 1.5 -- Merge the Training Cube with Mixture Model Cube
cube_merge_lsmm_train <- sits_merge(mm_cube, cube)

# ============================================================
# 2. Load and Explore Sample Data
# ============================================================

# 2.1 -- Load samples file (reference samples)
samples_sf  <- sf::st_read(file.path(sample_path, "samples_4_tiles-with-df-mask.gpkg"))

# 2.2 -- Extract Time Series from samples_sf and calculate the process duration
sits_get_data_start <- Sys.time()
samples <- sits_get_data(
  cube        = cube_merge_lsmm_train,
  samples     = samples_sf,
  label       = "label",
  n_sam_pol   = 16,       # number of pixels per segment
  multicores  = 8,       # parallel processing
  progress    = TRUE
)
sits_get_data_end <- Sys.time()
process_duration_sits_get_data <- round(sits_get_data_end-sits_get_data_start,2)
process_duration_sits_get_data

# 2.3.1 -- Visualize the temporal patterns of all features
plot(sits_patterns(samples))

# 2.3.2 -- Visualize the temporal patterns of specific features in a specific period
samples |> 
  sits_select(bands = c("SOIL","VEG","WATER"), start_date = '2023-08-01', end_date = '2025-07-28') |> 
  sits_patterns() |> 
  plot()

# 2.4 -- Save the samples Time Series to a R file
saveRDS(samples, paste0(rds_path, "/time_series/", process_version, tiles_train,  var, "_samples", ".rds"))

# ============================================================
# 3. 
# ============================================================

# 3.1 -- Clustering original Time Series Samples using SOM and calculate the process duration
# First, run with a 2x2 grid, then change to one of the values within the interval indicated by SITS and run again
sits_som_map_start <- Sys.time()
som_cluster <- sits_som_map(
  samples,
  grid_xdim = 33,
  grid_ydim = 33,
  alpha = 1.0,
  distance = "dtw",
  rlen = 20
)
sits_som_map_end <- Sys.time()
process_duration_sits_som_map <- round(sits_som_map_end-sits_som_map_start,2)
process_duration_sits_som_map

# 3.1.1 -- Plot the SOM map
plot(som_cluster)

ggsave(
  filename = paste0(process_version, tiles_train, var, "_", "som_eval.png"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# 3.1.2 -- Save the samples Time Series to a R file
saveRDS(som_cluster, paste0(rds_path, "/som/", process_version, "SOM1_33x33_", tiles_train, var, ".rds"))

# 3.1.3 -- Produce a tibble with a summary of the mixed labels:
som_eval <- sits_som_evaluate_cluster(som_cluster)

# 3.1.4 -- Plot the result of summary of the mixed labels
plot(som_eval) +
  labs(title = 'Confusão por cluster') +
  xlab("Porcentagem de Mistura") +
  ylab(NULL) +
  scale_fill_manual(values = my_colours, name = "Legenda") +
  theme(legend.position = "right")

ggsave(
  filename = paste0(process_version, tiles_train,  var, "_", "confusao_cluster.png"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# 3.2 -- Evaluates the quality of the samples based on the results of the SOM map
all_samples <- sits_som_clean_samples(
  som_map = som_cluster, 
  prior_threshold = 0.80,
  posterior_threshold = 0.70,
  keep = c("clean", "analyze", "remove"))

plot(all_samples)

ggsave(
  filename = paste0(process_version, tiles_train,  var, "_", "all_samples.png"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# 3.3 -- Filter samples according to evaluation 'clean' or to a specific class with low samples quantity 
clean_samples <- all_samples %>% filter(eval == "clean" | label == "DESMAT_DEGRAD_FOGO")

# 3.3.1 -- Save the new_samples Time Series to a R file
saveRDS(clean_samples, paste0(rds_path, "/time_series/", process_version, tiles_train,  var, "_clean_samples", ".rds"))

# 3.4 -- Clustering new Time Series Samples and calculate the process duration
# First, run with a 2x2 grid, then change to one of the values within the interval indicated by SITS and run again
sits_som_map_start2 <- Sys.time()
som_cluster_clean <- sits_som_map(data = clean_samples,
                                  grid_xdim = 32, 
                                  grid_ydim = 32,
                                  alpha = 1.0,
                                  distance = "dtw",
                                  rlen = 20
                                  # som_radius = 2,
                                  # mode = "online" # only for windows' PCs
                                  )
sits_som_map_end2 <- Sys.time()
process_duration_sits_som_map2 <- round(sits_som_map_end2-sits_som_map_start2,2)
process_duration_sits_som_map2

# 3.4.1 -- Plot the SOM map
plot(som_cluster_clean)

ggsave(
  filename = paste0(process_version, tiles_train,  var, "_", "som_eval_clean.png"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# 3.4.2 -- Save the new Time Series Samples to a R file
saveRDS(som_cluster_clean, paste0(rds_path, "/som/", process_version, "SOM2_32x32_", tiles_train,  var, ".rds"))

# 3.4.3 -- Produce a tibble with a summary of the mixed labels
som_eval_clean <- sits_som_evaluate_cluster(som_cluster_clean)

# 3.4.4 -- Plot the result
plot(som_eval_clean) +
  labs(title = 'Confusão por cluster') +
  xlab("Porcentagem de Mistura") +
  ylab(NULL) +
  scale_fill_manual(values = my_colours, name = "Legenda") +
  theme(legend.position = "right")

ggsave(
  filename = paste0(process_version, tiles_train,  var, "_", "confusao_cluster_clean.png"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# 3.5 -- Reduce imbalance between the classes and calculate the process duration
sits_reduce_imbalance_start <- Sys.time()
clean_samples_balanced <- sits_reduce_imbalance(
  samples = clean_samples,
  n_samples_over = 300,
  n_samples_under = 500
)
sits_reduce_imbalance_end <- Sys.time()
process_duration_sits_reduce_imbalance <- round(sits_reduce_imbalance_end-sits_reduce_imbalance_start,2)
process_duration_sits_reduce_imbalance

# 3.5.1 -- Removing columns that contain NA values
clean_samples_balanced <- clean_samples_balanced[, colSums(is.na(clean_samples_balanced)) == 0]

# 3.5.2 -- Save the new Time Series Samples Balanced to a R file
saveRDS(clean_samples_balanced, paste0(rds_path, "/time_series/", process_version, tiles_train,  var, "_clean_samples_balanced", ".rds"))

# 3.6 -- Clustering new Time Series Samples Balanced using SOM
# First, run with a 2x2 grid, then change to one of the values within the interval indicated by SITS and run again
sits_som_map_start3 <- Sys.time()
som_cluster_clean_balanced <- sits_som_map(clean_samples_balanced,
                                           grid_xdim = 17,
                                           grid_ydim = 17,
                                           alpha = 1.0,
                                           distance = "dtw",
                                           rlen = 20)
sits_som_map_end3 <- Sys.time()
process_duration_sits_som_map3 <- round(sits_som_map_end3-sits_som_map_start3,2)
process_duration_sits_som_map3

# 3.6.1 -- Plot the SOM map
plot(som_cluster_clean_balanced)

ggsave(
  filename = paste0(process_version, tiles_train,  var, "_", "som_eval_clean_balanced.png"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# 3.6.2 -- Produce a tibble with a summary of the mixed labels
som_eval_clean_balanced <- sits_som_evaluate_cluster(som_cluster_clean_balanced)

# 3.6.3 -- Plot the result
plot(som_eval_clean_balanced) +
  labs(title = 'Confusão por cluster') +
  xlab("Porcentagem de Mistura") +
  ylab(NULL) +
  scale_fill_manual(values = my_colours, name = "Legenda") +
  theme(legend.position = "right")

ggsave(
  filename = paste0(process_version, tiles_train,  var, "_", "confusao_cluster_clean_balanced.png"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# 3.6.4 -- Show the summary of the balanced time series sample data
summary(clean_samples_balanced)