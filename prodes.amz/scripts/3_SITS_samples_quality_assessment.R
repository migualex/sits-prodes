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
date_process <- format(Sys.Date(), "%Y-%m-%d_")
time_process <- format(Sys.time(), "%Hh%Mm", tz = "America/Sao_Paulo")
process_version <- paste0(date_process, time_process)

## III. Define the paths for files and folders needed in the processing
sample_path   <- "data/raw/samples/" #add the sample file to the path
rds_path      <- "data/rds/"
mixture_path  <- "data/raw/mixture_model"
plots_path    <- "data/plots/"

# IV. Identifier to distinguish this model run from previous versions
var <- "with-df-mask"

## V. Define a list with preference colors for each class
my_colors <- c(
  "OOB"                       = "black",
  "AGUA"                      = "#2980B9",
  "DESMAT_ARVORE_REMANESCE"   = "#a19c0a",
  "DESMAT_CORTE_RASO"         = "#f39c12",
  "DESMAT_CORTE_RASO_DM"      = "#f39c12",
  "DESMAT_DEGRAD_FOGO"        = "#EC7063",
  "DESMAT_VEG"                = "#D8DA83",
  "DESMAT_VEG_DM"             = "#D8DA83",
  "FLO_DEGRAD"                = "#9da676",
  "FLO_DEGRAD_FOGO"           = "#e6b0aa",
  "FLORESTA"                  = "#1E8449",
  "NF"                        = "#C0D665",
  "ROCHA"                     = "#229C59",
  "WETLANDS"                  = "#A0B9C8" 
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
  start_date  = "2022-08-01",
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
samples_sf  <- sf::st_read(file.path(sample_path, "samples-4-tiles-with-df-mask.gpkg"))

# 2.2 -- Extract Time Series from samples_sf and calculate the process duration
sits_get_data_start <- Sys.time()
samples <- sits_get_data(
  cube        = cube_merge_lsmm_train,
  samples     = samples_sf,
  label       = "label",
  n_sam_pol   = 16,       # number of pixels per segment
  multicores  = 8,       # adapt to your computer CPU core availability
  progress    = TRUE
)
sits_get_data_end <- Sys.time()
sits_get_data_time <- as.numeric(sits_get_data_end - sits_get_data_start, units = "secs")
sprintf("SITS get data process duration (HH:MM): %02d:%02d", as.integer(sits_get_data_time / 3600), as.integer((sits_get_data_time %% 3600) / 60))

# 2.3.1 -- Visualize the temporal patterns of all features
plot(sits_patterns(samples))

# 2.3.2 -- Visualize the temporal patterns of specific features in a specific period
samples |> 
  sits_select(bands = c("SOIL","VEG","WATER"), start_date = '2022-08-01', end_date = '2025-07-28') |> 
  sits_patterns() |> 
  plot()

# 2.4 -- Save the samples Time Series to a R file
saveRDS(samples, 
        paste0(rds_path,"time_series/", "samples-entire_", length(cube$tile),"-tiles-", tiles_train, "_", no.years,"-period-",cube_dates[1],"_",cube_dates[length(cube_dates)], "_", var, "_", process_version, ".rds"))


# ============================================================
# 3. Split samples in training and validation sets
# ============================================================

# 3.1 -- Set seed for reproducibility
set.seed(88)

message(sprintf("Total samples: %d", nrow(samples)))

# 3.2 -- Split training set (80%)
samples_train <- samples |> 
  group_by(label) |> 
  slice_sample(prop = 0.8) |> 
  ungroup()
message(sprintf("Training samples: %d (%.1f%%)", nrow(samples_train), 100 * nrow(samples_train) / nrow(samples)))

# 3.3 -- Creates the validation set by selecting samples NOT included in the training set (remaining 20%)
samples_val <- anti_join(samples, samples_train, by = names(samples))
message(sprintf("Validation samples: %d (%.1f%%)", nrow(samples_val), 100 * nrow(samples_val) / nrow(samples)))

# 3.4 -- Restore 'sits' class after dplyr operations
class(samples_train) <- c("sits", class(samples_train))
class(samples_val) <- c("sits", class(samples_val))

# 3.5 -- Check class distribution in training set
message("\n=== Training Set - Class Distribution ===")
print(table(samples_train$label))

# 3.6 -- Check class distribution in validation set
message("\n=== Validation Set - Class Distribution ===")
print(table(samples_val$label))

# 3.7 -- Visualize class distribution comparison
class_dist <- tibble(
  Class = names(table(samples$label)),
  Total = as.numeric(table(samples$label)),
  Training = as.numeric(table(samples_train$label)),
  Validation = as.numeric(table(samples_val$label))
) |>
  tidyr::pivot_longer(cols = c(Training, Validation), 
                      names_to = "Dataset", 
                      values_to = "Count")

p_split <- ggplot(class_dist, aes(x = Class, y = Count, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Training" = "#2ecc71", "Validation" = "#e74c3c")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Sample Distribution: Training vs Validation",
       x = "Class", y = "Number of Samples")

print(p_split)

# 3.8 -- Save the plot
ggsave(
  filename = paste0(plots_path, "train_val_split_", length(cube$tile),"-tiles-", tiles_train, "_", no.years,"-period-",cube_dates[1],"_",cube_dates[length(cube_dates)], "_", var, "_", process_version, ".png"),
  plot = p_split,
  width = 12,
  height = 6,
  dpi = 300
)

# 3.9 -- Save training and validation samples to RDS files
saveRDS(samples_train, 
        paste0(rds_path,"time_series/", "samples-train_", length(cube$tile),"-tiles-", tiles_train, "_", no.years,"-period-",cube_dates[1],"_",cube_dates[length(cube_dates)], "_", var, "_", process_version, ".rds"))

saveRDS(samples_val, 
        paste0(rds_path,"time_series/", "samples-val_", length(cube$tile),"-tiles-", tiles_train, "_", no.years,"-period-",cube_dates[1],"_",cube_dates[length(cube_dates)], "_", var, "_", process_version, ".rds"))


# ============================================================
# 4. Analyse quality, filter and balance
# ============================================================

# 4.1 -- Clustering original Time Series Samples using SOM and calculate the process duration
# First, run with a 2x2 grid, then change to one of the values within the interval indicated by SITS and run again
sits_som_map_start <- Sys.time()
som_cluster <- sits_som_map(
  samples_train,
  grid_xdim = 2,
  grid_ydim = 2,
  alpha = 1.0,
  distance = "dtw",
  rlen = 20
)
sits_som_map_end <- Sys.time()
sits_som_map_time <- as.numeric(sits_som_map_end - sits_som_map_start, units = "secs")
sprintf("SITS SOM map process duration (HH:MM): %02d:%02d", as.integer(sits_som_map_time / 3600), as.integer((sits_som_map_time %% 3600) / 60))

# 4.1.1 -- Plot the SOM map
plot(som_cluster)

ggsave(
  filename = paste0(process_version, "_", tiles_train, "_", var, "_som-eval.png"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# 4.1.2 -- Save the samples Time Series to a R file
saveRDS(som_cluster, paste0(
  rds_path, "som/",
  process_version,
  "_SOM1-", som_cluster$som$grid$xdim, "x", som_cluster$som$grid$ydim, "_",
  tiles_train, "_", var, ".rds")
)

# 4.1.3 -- Produce a tibble with a summary of the mixed labels:
som_eval <- sits_som_evaluate_cluster(som_cluster)

# 4.1.4 -- Plot the result of summary of the mixed labels
plot(som_eval) +
  labs(title = 'Confusion by cluster') +
  xlab("Percentage of Mixture") +
  ylab(NULL) +
  scale_fill_manual(values = my_colors, name = "Legend") +
  theme(legend.position = "right")

ggsave(
  filename = paste0(process_version, "_", tiles_train, "_", var, "_confusion-cluster.png"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# 4.2 -- Evaluates the quality of the samples based on the results of the SOM map
all_samples <- sits_som_clean_samples(
  som_map = som_cluster, 
  prior_threshold = 0.80,
  posterior_threshold = 0.70,
  keep = c("clean", "analyze", "remove"))

plot(all_samples)

ggsave(
  filename = paste0(process_version, "_", tiles_train, "_", var, "_all-samples.png"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# avaliar se a limpeza melhora ou não a qualidade da classificação. Testar com e sem a limpeza
# 4.3 -- Filter samples according to evaluation 'clean' or to a specific class with low samples quantity 
#clean_samples <- all_samples %>% filter(eval == "clean" | label == "DESMAT_DEGRAD_FOGO")
# 4.3.1 -- Save the new_samples Time Series to a R file
#saveRDS(clean_samples, paste0(rds_path, "time_series/", process_version, "_", tiles_train, "_", var, "_clean-samples.rds"))

# 4.4 -- Clustering new Time Series Samples and calculate the process duration
# First, run with a 2x2 grid, then change to one of the values within the interval indicated by SITS and run again
sits_som_map_start2 <- Sys.time()
som_cluster_clean <- sits_som_map(data = all_samples,
                                  grid_xdim = 2, 
                                  grid_ydim = 2,
                                  alpha = 1.0,
                                  distance = "dtw",
                                  rlen = 20
                                  # som_radius = 2,
                                  # mode = "online" # only for windows' PCs
                                  )
sits_som_map_end2 <- Sys.time()
sits_som_map_time2 <- as.numeric(sits_som_map_end2 - sits_som_map_start2, units = "secs")
sprintf("SITS SOM map 2 process duration (HH:MM): %02d:%02d", as.integer(sits_som_map_time2 / 3600), as.integer((sits_som_map_time2 %% 3600) / 60))

# 4.4.1 -- Plot the SOM map
plot(som_cluster_clean)

ggsave(
  filename = paste0(process_version, "_", tiles_train, "_", var, "_som-eval-clean.png"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# 4.4.2 -- Save the new Time Series Samples to a R file
saveRDS(som_cluster_clean, paste0(
  rds_path, "som/",
  process_version,
  "_SOM2-", som_cluster_clean$som$grid$xdim, "x", som_cluster_clean$som$grid$ydim, "_",
  tiles_train, "_", var, ".rds")
)

# 4.4.3 -- Produce a tibble with a summary of the mixed labels
som_eval_clean <- sits_som_evaluate_cluster(som_cluster_clean)

# 4.4.4 -- Plot the result
plot(som_eval_clean) +
  labs(title = 'Confusion by cluster') +
  xlab("Percentage of Mixture") +
  ylab(NULL) +
  scale_fill_manual(values = my_colors, name = "Legend") +
  theme(legend.position = "right")

ggsave(
  filename = paste0(process_version, "_", tiles_train, "_", var, "_confusion-cluster-clean.png"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# 4.5 -- Reduce imbalance between the classes and calculate the process duration
sits_reduce_imbalance_start <- Sys.time()
clean_samples_balanced <- sits_reduce_imbalance(
  samples = clean_samples,
  n_samples_over = 300,
  n_samples_under = 500
)
sits_reduce_imbalance_end <- Sys.time()
sits_reduce_imbalance_time <- as.numeric(sits_reduce_imbalance_end - sits_reduce_imbalance_start, units = "secs")
sprintf("SITS reduce imbalance process duration (HH:MM): %02d:%02d", as.integer(sits_reduce_imbalance_time / 3600), as.integer((sits_reduce_imbalance_time %% 3600) / 60))

# 4.5.1 -- Removing columns that contain NA values
clean_samples_balanced <- clean_samples_balanced[, colSums(is.na(clean_samples_balanced)) == 0]

# 4.5.2 -- Save the new Time Series Samples Balanced to a R file
saveRDS(clean_samples_balanced, paste0(rds_path, "time_series/", "samples-cleanned-&-balanced", "_", length(cube$tile),"-tiles-", tiles_train, "_", no.years,"-period-",cube_dates[1],"_",cube_dates[length(cube_dates)], "_", var, "_", process_version, ".rds"))

# 4.6 -- Clustering new Time Series Samples Balanced using SOM
# First, run with a 2x2 grid, then change to one of the values within the interval indicated by SITS and run again
sits_som_map_start3 <- Sys.time()
som_cluster_clean_balanced <- sits_som_map(clean_samples_balanced,
                                           grid_xdim = 2,
                                           grid_ydim = 2,
                                           alpha = 1.0,
                                           distance = "dtw",
                                           rlen = 20)
sits_som_map_end3 <- Sys.time()
sits_som_map_time3 <- as.numeric(sits_som_map_end3 - sits_som_map_start3, units = "secs")
sprintf("SITS SOM map 3 process duration (HH:MM): %02d:%02d", as.integer(sits_som_map_time3 / 3600), as.integer((sits_som_map_time3 %% 3600) / 60))

# 4.6.1 -- Plot the SOM map
plot(som_cluster_clean_balanced)

ggsave(
  filename = paste0(process_version, "_", tiles_train, "_", var, "_som-eval-clean-balanced.png"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# 4.6.2 -- Produce a tibble with a summary of the mixed labels
som_eval_clean_balanced <- sits_som_evaluate_cluster(som_cluster_clean_balanced)

# 4.6.3 -- Plot the result
plot(som_eval_clean_balanced) +
  labs(title = 'Confusion by cluster') +
  xlab("Percentage of Mixture") +
  ylab(NULL) +
  scale_fill_manual(values = my_colors, name = "Legend") +
  theme(legend.position = "right")

ggsave(
  filename = paste0(process_version, "_", tiles_train,  var, "_confusao-cluster-clean-balanced.png"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# 4.6.4 -- Show the summary of the balanced time series sample data
summary(clean_samples_balanced)