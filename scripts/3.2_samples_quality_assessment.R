# ============================================================
#  Samples quality assessment, filtering and balancing with SOM
# ============================================================

# ============================================================
# 1. Libraries, paths and some initial parameters
# ============================================================

# Step 1.1 -- Load Required Libraries
library(sits)
library(tibble)
library(dplyr)
library(ggplot2)

# Step 1.2 -- Define the date and time for the start of processing
date_process    <- format(Sys.Date(), "%Y-%m-%d_")
time_process    <- format(Sys.time(), "%Hh%Mm", tz = "America/Sao_Paulo")
process_version <- paste0(date_process, time_process)

# Step 1.3 -- Define the paths for files and folders needed in the processing
rds_path      <- "data/rds/"
plots_path    <- "data/plots/"
rds_filename  <- "samples_4-tiles-012014-012015-013014-013015_1y-period-2024-08-01_2025-07-28_prodes-amz_2026-02-24_10h30m.rds"

# Step 1.5 -- Identifier to distinguish this model run from previous versions
var <- "prodes-amz"

# Step 1.6 -- Define a list with preference colors for each class
my_colors <- c(
  "Water"                         = "#2980B9",
  "Wetland"                       = "#A0B9C8",
  "Forest"                        = "#1E8449",
  "Transition_Forest"             = "#E0DD22", 
  "Non_Forest_Natural_Vegetation" = "#C0D665",
  "Degradation"                   = "#9da676",
  "Degradation_Fire"              = "#e6b0aa",
  "Clear_Cut_Bare_Soil"           = "#f39c12",
  "Old_Clear_Cut_With_Vegetation" = "#B2B46D",
  "Clear_Cut_Burned_Area"         = "#CD6155",
  "Clear_Cut_With_Trees"          = "#a19c0a",
  "Clear_Cut_With_Vegetation"     = "#D8DA83",
  "Old_Clear_Cut_Bare_Soil"       = "#D39750"
)

# ATTENTION: Use the palette below if you are working in Non-Forest areas
my_colors <- c(
  "Hydrography_Lake"                                      = "#2980b9",
  "Hydrography_River"                                     = "#1f78b4",
  "Conversion_To_Agriculture"                             = "#f0b27a",
  "Conversion_To_Bare_Soil"                               = "#f39c12",
  "Previous_Conversion_To_Agriculture"                    = "#b08b57",
  "Previous_Conversion_To_Bare_Soil"                      = "#a0522d",
  "Fire_In_Non_Forest_Natural_Vegetation"                 = "#cd6155",
  "Non_Forest_Natural_Vegetation_Dry_High_Biomass"        = "#f6cc41",
  "Non_Forest_Natural_Vegetation_Dry_Low_Biomass"         = "#dbebd8",
  "Non_Forest_Natural_Vegetation_Post_Fire"               = "#e6b0aa",
  "Non_Forest_Natural_Vegetation_Wet"                     = "#a0b9c8",
  "Non_Forest_Natural_Vegetation_Forest_Transition"       = "#88cda2",
  "Non_Forest_Natural_Vegetation_Woodland"                = "#1e8449",
  "Non_Forest_Natural_Vegetation_Palm_Swamps"             = "#3ababa"
)

# ============================================================
# 3. Load and Explore Train Sample Data
# ============================================================

# Step 3.3 -- Load the samples Time Series from a R file
samples <- readRDS(file.path(rds_path, "time_series", rds_filename))

# Step 3.3 -- Create output directory per tile and period
tiles_train <- gsub(".*_(\\d{6}(-\\d{6})*)_.*", "\\1", rds_filename)

tile_period_dir <- file.path(plots_path, var)
dir.create(tile_period_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 4. Analyse quality (SOM - 1)
# ============================================================

# Step 4.1 -- Clustering original Time Series Samples using SOM and calculate the process duration
# First, run with a 2x2 grid, then change to one of the values within the interval indicated by SITS and run again
sits_som_map_start <- Sys.time()
som_cluster <- sits_som_map(
  samples,
  grid_xdim = 2,
  grid_ydim = 2,
  alpha = 1.0,
  distance = "dtw",
  rlen = 20
)
sits_som_map_end <- Sys.time()
sits_som_map_time <- as.numeric(sits_som_map_end - sits_som_map_start, units = "secs")
sprintf("SITS SOM map process duration (HH:MM): %02d:%02d", as.integer(sits_som_map_time / 3600), as.integer((sits_som_map_time %% 3600) / 60))

# Step 4.1.1 -- Plot the SOM map 1
plot(som_cluster)

# Step 4.1.2 -- Save SOM map plot
ggsave(
  filename = paste0(process_version, "_", tiles_train, "_", var, "_som-eval.png"),
  path = tile_period_dir,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# Step 4.2 -- Save the samples Time Series to a R file
saveRDS(som_cluster, paste0(
  rds_path, "som/",
  process_version,
  "_SOM1-", som_cluster$som$grid$xdim, "x", som_cluster$som$grid$ydim, "_",
  tiles_train, "_", var, ".rds")
)

# Step 4.3 -- Produce a tibble with a summary of the mixed labels:
som_eval <- sits_som_evaluate_cluster(som_cluster)

# Step 4.3.1 -- Plot the result of summary of the mixed labels
plot(som_eval) +
  labs(title = 'Confusion by cluster') +
  xlab("Percentage of Mixture") +
  ylab(NULL) +
  scale_fill_manual(values = my_colors, name = "Legend") +
  theme(legend.position = "right")

# Step 4.3.2 -- Save the plot of summary of the mixed labels
ggsave(
  filename = paste0(process_version, "_", tiles_train, "_", var, "_confusion-cluster.png"),
  path = tile_period_dir,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# ============================================================
# 5. Analyse quality, filter (SOM - 2)
# ============================================================

# Step 5.1 -- Evaluates the quality of the samples based on the results of the SOM map
all_samples <- sits_som_clean_samples(
  som_map = som_cluster, 
  prior_threshold = 0.80,
  posterior_threshold = 0.70,
  keep = c("clean", "analyze", "remove"))

# Step 5.1.1 -- Plot quality of the samples
plot(all_samples)

# Step 5.1.2 -- Save the plot
ggsave(
  filename = paste0(process_version, "_", tiles_train, "_", var, "_all-samples.png"),
  path = tile_period_dir,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# Step 5.3 -- Filter samples according to evaluation 'clean' or to a specific class with low samples quantity
clean_samples <- all_samples %>% filter(eval == "clean" | label == "DESMAT_DEGRAD_FOGO")

# Step 5.3.1 -- Save the new_samples Time Series to a R file
saveRDS(clean_samples, paste0(rds_path, "time_series/", process_version, "_", tiles_train, "_", var, "_clean-samples.rds"))

# Step 5.4 -- Clustering new Time Series Samples and calculate the process duration
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

# Step 5.4.1 -- Plot the SOM map 2
plot(som_cluster_clean)

# Step 5.4.2 -- Save the plot
ggsave(
  filename = paste0(process_version, "_", tiles_train, "_", var, "_som-eval-clean.png"),
  path = tile_period_dir,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# Step 5.4.3 -- Save the new Time Series Samples to a R file
saveRDS(som_cluster_clean, paste0(
  rds_path, "som/",
  process_version,
  "_SOM2-", som_cluster_clean$som$grid$xdim, "x", som_cluster_clean$som$grid$ydim, "_",
  tiles_train, "_", var, ".rds")
)

# Step 5.5 -- Produce a tibble with a summary of the mixed labels
som_eval_clean <- sits_som_evaluate_cluster(som_cluster_clean)

# Step 5.5.1 -- Plot the result
plot(som_eval_clean) +
  labs(title = 'Confusion by cluster') +
  xlab("Percentage of Mixture") +
  ylab(NULL) +
  scale_fill_manual(values = my_colors, name = "Legend") +
  theme(legend.position = "right")

# Step 5.5.2 -- Save the plot
ggsave(
  filename = paste0(process_version, "_", tiles_train, "_", var, "_confusion-cluster-clean.png"),
  path = tile_period_dir,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# ============================================================
# 6. Analyse quality, filter and balance (SOM - 3)
# ============================================================

# Step 6.1 -- Reduce imbalance between the classes and calculate the process duration
sits_reduce_imbalance_start <- Sys.time()
clean_samples_balanced <- sits_reduce_imbalance(
  samples = clean_samples,
  n_samples_over = 300,
  n_samples_under = 500
)
sits_reduce_imbalance_end <- Sys.time()
sits_reduce_imbalance_time <- as.numeric(sits_reduce_imbalance_end - sits_reduce_imbalance_start, units = "secs")
sprintf("SITS reduce imbalance process duration (HH:MM): %02d:%02d", as.integer(sits_reduce_imbalance_time / 3600), as.integer((sits_reduce_imbalance_time %% 3600) / 60))

# Step 6.2 -- Removing columns that contain NA values
clean_samples_balanced <- clean_samples_balanced[, colSums(is.na(clean_samples_balanced)) == 0]

# Step 6.2.1 -- Save the new Time Series Samples Balanced to a R file
saveRDS(clean_samples_balanced, paste0(rds_path, "time_series/", "samples-cleanned-&-balanced", "_", tiles_train, "_", var, "_", process_version, ".rds"))

# Step 6.3 -- Clustering new Time Series Samples Balanced using SOM
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

# Step 6.3.1 -- Plot the SOM map
plot(som_cluster_clean_balanced)

# Step 6.3.2 -- Save the plot
ggsave(
  filename = paste0(process_version, "_", tiles_train, "_", var, "_som-eval-clean-balanced.png"),
  path = tile_period_dir,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# Step 6.4 -- Produce a tibble with a summary of the mixed labels
som_eval_clean_balanced <- sits_som_evaluate_cluster(som_cluster_clean_balanced)

# Step 6.4.1-- Show the summary of the balanced time series sample data
summary(clean_samples_balanced)

# Step 6.4.2-- Plot the result
plot(som_eval_clean_balanced) +
  labs(title = 'Confusion by cluster') +
  xlab("Percentage of Mixture") +
  ylab(NULL) +
  scale_fill_manual(values = my_colors, name = "Legend") +
  theme(legend.position = "right")

# Step 6.4.3 -- Save the plot
ggsave(
  filename = paste0(process_version, "_", tiles_train,  var, "_confusao-cluster-clean-balanced.png"),
  path = tile_period_dir,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)
print("Samples analyzed successfully!")
