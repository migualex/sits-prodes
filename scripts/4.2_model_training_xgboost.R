# ============================================================
# Train a XGBoost machine learning model
# ============================================================

# Load required libraries
library(sits)
library(ggplot2)

# Define the parameters: These are user-defined variables
time_series_name  <- "TS-tiles_012014-012015-013014-013015_1y_2024-08-01_2025-07-31_all-samples-new-pol-avg-false_2026-04-22_10h46m.rds"
start_date        <- "2024-08-01"
end_date          <- "2025-07-31"
tiles             <- c("012014","012015","013014","013015")

# Function to read class names and their colors::IMPORTANT
read_class_config <- function(config_file = "class_config.txt") {
  
  if (!file.exists(config_file)) {
    stop(paste("Configuration file not found:", config_file))
  }
  
  lines <- readLines(config_file, encoding = "UTF-8", warn = FALSE)
  
  # Remove empty lines and comments
  lines <- trimws(lines)
  lines <- lines[nchar(lines) > 0 & !startsWith(lines, "#")]
  
  # Identify sections and populate lists
  current_section  <- NULL
  class_trans_list <- list()
  colors_list      <- list()
  
  for (line in lines) {
    if (startsWith(line, "[") && endsWith(line, "]")) {
      current_section <- gsub("\\[|\\]", "", line)
      next
    }
    
    if (!is.null(current_section) && grepl("=", line)) {
      parts <- strsplit(line, "=", fixed = TRUE)[[1]]
      key   <- trimws(parts[1])
      value <- trimws(paste(parts[-1], collapse = "=")) # preserves '=' in hex codes
      
      if (current_section == "CLASS_TRANSLATION") {
        class_trans_list[[key]] <- value
      } else if (current_section == "COLORS") {
        colors_list[[key]] <- value
      }
    }
  }
  
  class_translation <- unlist(class_trans_list)
  my_colors         <- unlist(colors_list)
  
  message(sprintf("Config loaded: %d class translations | %d colors",
                  length(class_translation), length(my_colors)))
  
  return(list(
    class_translation = class_translation,
    my_colors         = my_colors
  ))
}

# Date and time of the start of processing
date_process    <- format(Sys.Date(), "%Y-%m-%d_")
time_process    <- format(Sys.time(), "%Hh%Mm", tz = "America/Sao_Paulo")
process_version <- paste0(date_process, time_process)

# File and folder paths
time_series_path  <- file.path("data/rds/time_series/", time_series_name)
rds_path          <- "data/rds/"
plots_path        <- "data/plots/"
config_dir        <- ".."

# Identifier to distinguish this model run from previous versions
var <- stringr::str_split_i(time_series_name, "_", 6)

# ============================================================
# 1. Define and Load Data Cubes
# ============================================================

# Step 1.1 -- Create a training cube from a collection
cube <- sits_cube(
  source      = "BDC",
  collection  = "SENTINEL-2-16D",
  bands       = c('B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12', 'NDVI', 'NBR', 'EVI', 'CLOUD'),
  tiles       = tiles,
  start_date  = start_date,
  end_date    = end_date,
  progress    = TRUE)

# Step 1.2 -- Calculate the number of years in the training cube
no.years <- paste0(floor(lubridate::year(end_date) - lubridate::year(start_date)), "y")
tiles_train <- paste(sort(tiles), collapse = "-")
no.cubes <- paste0(length(cube$tile), "t")

# ============================================================
# 2. Cross-validation of training data
# ============================================================

# Step 2.1 -- Reading training samples
train_samples <- readRDS(time_series_path)

# Step 2.2 -- Load color palette from external config file
config     <- read_class_config(file.path(config_dir, "class_config.txt"))
my_colors  <- config$my_colors
my_colors  <- my_colors[names(my_colors) %in% unique(train_samples$label)]

# ============================================================
# 3. Training and saving model
# ============================================================

# Step 3.1 -- Set random seed to ensure reproducibility
set.seed(88)

# Record start time of the hyperparameter tuning process
sits_xgb_model_fine_tuning_start <- Sys.time()

# Perform hyperparameter tuning for the XGBoost model
tunned_model <- sits_tuning(
  train_samples,
  validation_split = 0.3,
  ml_method = sits_xgboost(),
  params = sits_tuning_hparams(
    learning_rate = loguniform(0.01, 0.3),        # learning rate (log-scale search)
    min_split_loss = uniform(1,10),                # gamma: minimum loss reduction to split
    max_depth = choice(6,7,8,9,10,11,12,13,14,15), # maximum tree depth
    min_child_weight = choice(5, 10, 15, 20, 25, 30), # minimum sum of instance weight in a leaf
    subsample = 1,                      # fraction of samples used per tree (1 = no subsampling)
    nrounds = choice(100, 200, 300, 400, 500, 600, 700), # number of boosting rounds
    nthread = 28                        # number of threads (parallel processing)
  ),
  trials = 30,       # number of parameter combinations tested
  multicores = 28,   # parallel tuning
  progress = TRUE    # display progress
)

# Record end time and compute total duration
sits_xgb_model_fine_tuning_end <- Sys.time()
sits_xgb_model_fine_tuning_time <- as.numeric(
  sits_xgb_model_fine_tuning_end - sits_xgb_model_fine_tuning_start,
  units = "secs"
)

# Print duration in HH:MM format
sprintf(
  "SITS XGBoost model fine tuning process duration (HH:MM): %02d:%02d",
  as.integer(sits_xgb_model_fine_tuning_time / 3600),
  as.integer((sits_xgb_model_fine_tuning_time %% 3600) / 60)
)

# Flatten list-columns of hyperparameters into atomic vectors
tunned_model <- tunned_model |>
  dplyr::mutate(
    dplyr::across(learning_rate:verbose, ~ unlist(.x))
)

# Step 3.2 -- Build confusion matrix as a tibble
matriz_conf_xbb_model <- tibble::tibble(
  as.data.frame(tunned_model$acc[[1]]$table)
)

# Step 3.3 -- Plot confusion matrix
ggplot(matriz_conf_xbb_model, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 4) +
  scale_fill_distiller(
    palette = "Blues",
    direction = 1,
    name = "Cases"
  ) +
  labs(
    x = "Reference",
    y = "Predicted", 
    title = "Confusion Matrix"
  ) +
  scale_y_discrete(
    limits = rev,
    labels = gsub("_", " ", config$class_translation)
  ) +
  scale_x_discrete(
    labels = gsub("_", " ", config$class_translation)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, size = 8), 
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )

# Step 3.4 -- Train final model using best hyperparameters
xgb_model <- sits_train(
  samples = train_samples,
  ml_method = sits_xgboost(
    learning_rate = tunned_model$learning_rate[1],
    min_split_los = tunned_model$min_split_loss[1], # NOTE: check typo (min_split_los vs min_split_loss)
    max_depth = tunned_model$max_depth[1],
    min_child_weight = tunned_model$min_child_weight[1],
    subsample = 1,
    nrounds = tunned_model$nrounds[1],
    nthread = 28,
    verbose = TRUE
  )
)

# Step 3.4.1 -- Plot feature importance of the trained model
plot(xgb_model)

# Step 3.4.2 -- Save the generated plot
ggsave(
  filename = paste0(
    process_version, "_", tiles_train, "_", no.years, var,
    "_minimal_tree_depth.png"
  ),
  path = tile_period_dir,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350
)

# Step 3.5 -- Save trained model to disk (RDS format)
saveRDS(
  xgb_model,
  paste0(
    rds_path, "model/xgboost/", "XGB-model_",
    length(cube$tile), "-tiles-", tiles_train, "_",
    no.years, "-period-",
    cube_dates[1], "_", cube_dates[length(cube_dates)],
    "_", var, "_", process_version, ".rds"
  )
)
print("Model trained successfully!")
