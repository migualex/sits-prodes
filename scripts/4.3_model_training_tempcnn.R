# ============================================================
# Train a Temporal Convolution Neural Network deep learning model
# ============================================================

# Load required libraries
library(sits)
library(ggplot2)
library(torch)
library(luz)

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

# Step 2.3 -- Using k-fold validation
sits_kfold_validate_start <- Sys.time()
rfor_validate <- sits_kfold_validate(
  samples = train_samples,
  folds = 5, # how many times to split the data (default = 5)
  ml_method = sits_rfor(),
  multicores = 28,
  progress = TRUE) # adapt to your computer CPU core availability
sits_kfold_validate_end <- Sys.time()
sits_kfold_validate_time <- as.numeric(sits_kfold_validate_end - sits_kfold_validate_start, units = "secs")
sprintf("SITS kfold_validate process duration (HH:MM): %02d:%02d", as.integer(sits_kfold_validate_time / 3600), as.integer((sits_kfold_validate_time %% 3600) / 60))

# Step 2.3.1 -- Plot the confusion matrix
plot(rfor_validate, type = "confusion_matrix")

# Step 2.3.2 -- Plot the metrics by class
plot(rfor_validate, type = "metrics")

# ============================================================
# 3. Training and saving model
# ============================================================

# Step 3.1 -- Set random seed to ensure reproducibility
set.seed(88)

# Step 3.2 -- Train TempCNN model using training samples
tempcnn_model <- sits_train(
  samples = train_samples,
  ml_method = sits_tempcnn(
    
    # Convolutional layers configuration
    cnn_layers = c(64, 64, 64),              # number of filters per layer
    cnn_kernels = c(3, 3, 3),                # kernel size for each convolution
    cnn_dropout_rates = c(0.2, 0.2, 0.2),    # dropout rate to reduce overfitting
    
    # Fully connected (dense) layer configuration
    dense_layer_nodes = 256,                 # number of neurons in dense layer
    dense_layer_dropout_rate = 0.5,          # dropout for dense layer
    
    # Training configuration
    epochs = 150,                            # number of training epochs
    batch_size = 1024,                       # batch size for training
    validation_split = 0.2,                  # proportion of data used for validation
    
    # Optimizer configuration (AdamW)
    optimizer = torch::optim_adamw,
    opt_hparams = list(
      lr = 5e-04,                            # learning rate
      eps = 1e-08,                           # numerical stability term
      weight_decay = 1e-06                   # L2 regularization
    ),
    
    # Learning rate scheduling
    lr_decay_epochs = 1,                     # frequency of decay (in epochs)
    lr_decay_rate = 0.95,                    # decay factor
    
    # Early stopping configuration
    patience = 20,                           # epochs to wait without improvement
    min_delta = 0.01,                        # minimum improvement threshold
    
    verbose = FALSE                          # disable training logs
  )
)

# Step 3.2.1 -- Plot model diagnostics (e.g., training history or feature relevance)
plot(tempcnn_model)

# Step 3.3 -- Export trained model to a more accessible structure
tempcnn_model2 <- sits_model_export(tempcnn_model)

# Step 3.3.1 -- Open PNG device to save diagnostic plots
png(
  filename = file.path(
    tile_period_dir,
    paste0(
      process_version, "_", tiles_train, "_",
      no.years, var, "_oob_ntree_mde.png"
    )
  ),
  width = 3529,
  height = 1578,
  res = 350
)

# Step 3.3.2 -- Plot Out-Of-Bag (OOB) error over number of trees
matplot(
  tempcnn_model2$err.rate, 
  type = "l", lty = 1, lwd = 2,
  col = my_colors,           
  main = "Out-of-Bag Error vs Number of Trees",
  xlab = "Number of Trees (ntree)", 
  ylab = "Out-of-Bag Error"
)

# Step 3.3.3 -- Add legend to the plot
legend(
  "topright", 
  legend = names(my_colors), 
  col = my_colors, 
  lty = 1,      
  cex = 1,    
  bty = "n"
)

# Close PNG device and save file
dev.off()

# Step 3.4 -- Save trained model to disk (RDS format)
saveRDS(
  tempcnn_model,
  paste0(
    rds_path, "model/tempcnn/", "TCNN-model_",
    length(cube$tile), "-tiles-", tiles_train, "_",
    no.years, "-period-",
    cube_dates[1], "_", cube_dates[length(cube_dates)],
    "_", var, "_", process_version, ".rds"
  )
)
print("Model trained successfully!")
