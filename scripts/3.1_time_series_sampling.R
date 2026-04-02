# ============================================================
#  Time series extraction from the training samples
# ============================================================

# ============================================================
# 1. Libraries, paths and some initial parameters
# ============================================================

# Step 1.1 -- Load Required Libraries
library(sits)
library(tibble)
library(dplyr)
library(ggplot2)

# Step 1.2 -- Function to read class names and their colors::IMPORTANT
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

# Step 1.3 -- Define the date and time for the start of processing
date_process    <- format(Sys.Date(), "%Y-%m-%d_")
time_process    <- format(Sys.time(), "%Hh%Mm", tz = "America/Sao_Paulo")
process_version <- paste0(date_process, time_process)

# Step 1.4 -- Define the paths for files and folders needed in the processing
sample_path   <- "data/raw/samples" #add the sample file to the path
rds_path      <- "data/rds/"
mixture_path  <- "data/raw/mixture_model"
plots_path    <- "data/plots/"
config_dir    <- "../scripts"

# Step 1.5 -- Define time range
start_date   <- "2024-08-01"
end_date     <- "2025-07-31"
tiles        <- c("012014","012015","013014","013015")

# Step 1.6 -- Identifier to distinguish this model run from previous versions
var <- "prodes-amz"

# ============================================================
# 2. Define and Load Data Cubes
# ============================================================

# Step 2.1 -- Create a training cube from a collection
cube <- sits_cube(
  source      = "BDC",
  collection  = "SENTINEL-2-16D",
  bands       = c('B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12', 'NDVI', 'NBR', 'EVI', 'CLOUD'),
  tiles       = tiles,
  start_date  = start_date,
  end_date    = end_date,
  progress    = TRUE)

# Step 2.2 -- Calculate the number of years in the training cube
cube_dates <- sits_timeline(cube)
no.years <- paste0(floor(lubridate::year(end_date) - lubridate::year(start_date)), "y")

# Step 2.3 -- Concatenates all the names of the training tiles into a single string separated by '-'
tiles_train <- paste(cube$tile, collapse = "-")

# Step 2.4 -- Retrieve Mixture Model Cube from a predefined repository
mm_cube <- sits_cube(
  source      = "BDC",
  collection  = "SENTINEL-2-16D",
  bands       = c("SOIL", "VEG", "WATER"),
  tiles       = tiles,
  data_dir    = mixture_path,
  start_date  = start_date,
  end_date    = end_date,
  progress    = TRUE)

# Step 2.5 -- Merge the Training Cube with Mixture Model Cube
cube_merge_lsmm_train <- sits_merge(mm_cube, cube)

# Step 2.6 -- Create output directory per tile and period
tiles_id <- paste(sort(unique(cube_merge_lsmm_train$tile)), collapse = "_")

tile_period_dir <- file.path(plots_path, var)

dir.create(tile_period_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 3. Load and Explore Train Sample Data
# ============================================================

# Step 3.1 -- Read training samples (rewrite the name of your samples file)
sampling_date   <- "2026-02-24"   # date of the sampling file (YYYY-MM-DD)
tiles_str       <- paste(sort(tiles), collapse = "-")
samples_name    <- paste("sampling-training", tiles_str, var, sampling_date, sep = "-")
samples_train   <- sf::st_read(file.path(sample_path, paste0(samples_name, ".gpkg")))

# Step 3.2 -- Load class translation from external config file
config     <- read_class_config(file.path(config_dir, "class_config.txt"))
class_translation <- config$class_translation

# Step 3.2.1 -- Apply translation: keep the original label if no translation is found
samples_train$label <- ifelse(
  samples_train$label %in% names(class_translation),
  class_translation[samples_train$label],
  samples_train$label
)

print(table(samples_train$label))

# Step 3.3 -- Extract Time Series from samples_train and calculate the process duration
sits_get_data_start <- Sys.time()
samples <- sits_get_data(
  cube        = cube_merge_lsmm_train,
  samples     = samples_train,
  n_sam_pol   = 16,
  pol_avg     = FALSE,
  label       = "label",
  multicores  = 28,       # adapt to your computer CPU core availability
  progress    = TRUE)
sits_get_data_end <- Sys.time()
sits_get_data_time <- as.numeric(sits_get_data_end - sits_get_data_start, units = "secs")
sprintf("SITS get data process duration (HH:MM): %02d:%02d", as.integer(sits_get_data_time / 3600), as.integer((sits_get_data_time %% 3600) / 60))

# Step 3.3.1 -- Visualize the temporal patterns of all features
plot(sits_patterns(samples))

# Step 3.3.2 -- Visualize the temporal patterns of specific features in a specific period
samples |> 
  sits_select(bands = c("NDVI","B04","B08","B11"), start_date = '2024-08-12', end_date = '2025-07-28') |> 
  sits_patterns() |> 
  plot()

# Step 3.4 -- Save the samples Time Series to a R file
saveRDS(samples, 
        paste0(rds_path,"time_series/", "samples_", length(cube$tile),"-tiles-", tiles_train, "_", no.years,"-period-",cube_dates[1],"_",cube_dates[length(cube_dates)], "_", var, "_", process_version, ".rds"))
print("Time series extracted successfully!")
