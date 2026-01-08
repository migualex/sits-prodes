# ============================================================
#  Random Forest model training
# ============================================================

## I. Load Required Libraries
library(sits)
library(sf)
library(tibble)
library(ggplot2)
library(lubridate)

## II. Define the date and time for the start of processing
date_process <- format(Sys.Date(), "/%Y_%m_%d_")
time_process <- format(Sys.time(), "%Hh%Mm_", tz = "America/Sao_Paulo")
process_version <- paste0(date_process, time_process)

## III. Define the paths for files and folders needed in the processing
sample_path   <- "~/sits_amz/samples/AMOSTRAS.gpkg"
rds_path      <- "~/sits_amz/rds/2025_12_15"
mixture_path  <- "~/sits_amz/mixture/all_tiles_1y"
plots_path    <- ""

## IV. Define a list with preference colours for each class
my_colours <- c(
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
  start_date  = "2024-08-01",
  end_date    = "2025-07-31",
  progress    = TRUE
)

datas <- sits_timeline(cube)
qtd_anos <- paste0(floor(interval(datas[1], datas[length(datas)]) / years(1)), "y")

# Concatenates all the names of the training tiles into a single string separated by '_'
tiles_train <- paste(cube$tile, collapse = "_")

# Step 1.2 -- Retrieve Mixture Model Cube from a predefined repository
mm_cube <- sits_cube(
  source = "BDC",
  tiles = c('012014',  '012015', '013014', '013015'),
  collection = "SENTINEL-2-16D",
  bands = c("SOIL", "VEG", "WATER"),
  data_dir = mixture_path,
  progress = TRUE
)

# Merge the Training Cube with Mixture Model Cube
cube_merge_lsmm_train <- sits_merge(mm_cube, cube)

# ============================================================
# 2. Load and Explore Sample Data
# ============================================================

# 2.1 -- Load samples file (reference samples)
samples_sf  <- sf::st_read(sample_path)


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
  sits_select(bands = c("SOIL","VEG","WATER"), start_date = '2024-08-01', end_date = '2025-07-28') |> 
   sits_patterns() |> 
   plot()


# 2.4 -- Save the samples Time Series to a R file
saveRDS(samples,paste0(rds_path, process_version, "TS_", tiles_train,".rds"))


# ============================================================
# 3. Samples quality assessment, filtering and balancing
# ============================================================

# 3.1 -- Clustering original Time Series Samples using SOM and calculate the process duration
sits_som_map_start <- Sys.time()
som_cluster <- sits_som_map(
  samples,
  grid_xdim = 15,
  grid_ydim = 15,
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
  filename = paste0(process_version, tiles_train, "", "som_eval.png"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# 3.1.2 -- Save the samples Time Series to a R file
saveRDS(som_cluster,paste0(rds_path, process_version, "SOM1_15x15_", tiles_train,".rds"))

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
  filename = paste0(process_version, tiles_train, "", "confusao_cluster.png"),
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
  prior_threshold = 0.7,
  posterior_threshold = 0.6,
  keep = c("clean", "analyze", "remove"))

plot(all_samples)

ggsave(
  filename = paste0(process_version, tiles_train, "", "all_samples.png"),
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
saveRDS(clean_samples,paste0(rds_path, process_version,"clean_samples_", tiles_train,".rds"))


# 3.4 -- Clustering new Time Series Samples and calculate the process duration
sits_som_map_start2 <- Sys.time()
som_cluster_clean <- sits_som_map(data = clean_samples,
                                grid_xdim = 15,
                                grid_ydim = 15,
                                alpha = 1.0,
                                distance = "dtw",
                                rlen = 20)
                                # som_radius = 2,
                                # mode = "online" apenas em PC com windows)
sits_som_map_end2 <- Sys.time()
process_duration_sits_som_map2 <- round(sits_som_map_end2-sits_som_map_start2,2)
process_duration_sits_som_map2

# 3.4.1 -- Plot the SOM map
plot(som_cluster_clean)

ggsave(
  filename = paste0(process_version, tiles_train, "", "som_eval_clean.png"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# 3.4.2 -- Save the new Time Series Samples to a R file
saveRDS(som_cluster_clean,paste0(rds_path, process_version, "SOM2_15x15_", tiles_train,".rds"))

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
  filename = paste0(process_version, tiles_train, "", "confusao_cluster_clean.png"),
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
  n_samples_under = 1000
)
sits_reduce_imbalance_end <- Sys.time()
process_duration_sits_reduce_imbalance <- round(sits_reduce_imbalance_end-sits_reduce_imbalance_start,2)
process_duration_sits_reduce_imbalance

# 3.5.1 -- Removing columns that contain NA values
clean_samples_balanced <- clean_samples_balanced[, colSums(is.na(clean_samples_balanced)) == 0]

# 3.5.2 -- Save the new Time Series Samples Balanced to a R file
saveRDS(clean_samples_balanced,paste0(rds_path, process_version,"clean_samples_balanced", tiles_train,".rds"))


# 3.6 -- Clustering new Time Series Samples Balanced using SOM
sits_som_map_start3 <- Sys.time()
som_cluster_clean_balanced <- sits_som_map(clean_samples_balanced,
                                          grid_xdim = 15,
                                          grid_ydim = 15,
                                          alpha = 1.0,
                                          distance = "dtw",
                                          rlen = 20)
sits_som_map_end3 <- Sys.time()
process_duration_sits_som_map3 <- round(sits_som_map_end3-sits_som_map_start3,2)
process_duration_sits_som_map3

# 3.6.1 -- Plot the SOM map
plot(som_cluster_clean_balanced)

ggsave(
  filename = paste0(process_version, tiles_train, "", "som_eval_clean_balanced.png"),
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
  filename = paste0(process_version, tiles_train, "", "confusao_cluster_clean_balanced.png"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# 3.6.4 -- Show the summary of the balanced time series sample data
summary(clean_samples_balanced)


# ============================================================
# 4. Train Classification Model
# ============================================================

# Step 4.1 -- Set a seed of random number generator (RNG) for reproducibility
set.seed(88)

# Step 4.2 -- Train the Random Forest model
rf_model <- sits_train(
   samples   = clean_samples_balanced,
   ml_method = sits_rfor()
 )

# Step 4.3 -- Plot the most important variables of the model
plot(rf_model)

var <- stringr::str_extract(basename(sample_path), "_(with|no)_df_mask")

# Step 4.4 -- Save the model to a R file
saveRDS(rf_model,paste0(rds_path, process_version,"RF_", qtd_anos, "_", tiles_train, var,".rds"))

Print("Model has been trained!")