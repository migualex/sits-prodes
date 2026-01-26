# ============================================================
# Train Machine Learning Model
# ============================================================

## I. Load Required Libraries
library(sits)
library(ggplot2)

## II. Define the date and time for the start of processing
date_process <- format(Sys.Date(), "%Y-%m-%d_")
time_process <- format(Sys.time(), "%Hh%Mm", tz = "America/Sao_Paulo")
process_version <- paste0(date_process, time_process)

## III. Define the paths for files and folders needed in the processing
rds_path    <- "data/rds/"
plots_path  <- "data/plots/"

# IV. Identifier to distinguish this model run from previous versions | Keep it the same as the samples' identifier
var <- "with-df-mask-with-all-samples"

## V. Define a list with preference colours for each class
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


# ============================================================
# 2. Cross-validation of training data
# ============================================================

# Step 2.1 -- Reading training samples
train_samples <- readRDS("~/sits-prodes/prodes.amz/data/rds/time_series/2026-01-21_10h10m_012015-012014-013015-013014_with-df-mask_clean-samples.rds")

# Step 2.2 -- Using k-fold validation
sits_kfold_validate_start <- Sys.time()
rfor_validate <- sits_kfold_validate(
  samples = train_samples,
  folds = 5, # how many times to split the data (default = 5)
  ml_method = sits_rfor(),
  multicores = 12,
  progress = TRUE) # adapt to your computer CPU core availability
sits_kfold_validate_end <- Sys.time()
sits_kfold_validate_time <- as.numeric(sits_kfold_validate_end - sits_kfold_validate_start, units = "secs")
sprintf("SITS kfold_validate process duration (HH:MM): %02d:%02d", as.integer(sits_kfold_validate_time / 3600), as.integer((sits_kfold_validate_time %% 3600) / 60))

# Step 2.2.1 -- Plot the confusion matrix
plot(rfor_validate, type = "confusion_matrix")

# Step 2.2.2 -- Plot the metrics by class
plot(rfor_validate, type = "metrics")


# ============================================================
# 3. Training and saving model
# ============================================================

# Step 3.1 -- Set a seed of random number generator (RNG) for reproducibility
set.seed(88)

# Step 3.2 -- Train the model
rf_model <- sits_train(
   samples   = train_samples,
   ml_method = sits_rfor(num_trees = 100)
 )

# Step 3.3.1 -- Plot the most important variables of the model
plot(rf_model)

ggsave(
  filename = paste0(process_version, "_", tiles_train,"_", no.years, var, "_minimal_tree_depth.png"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)


# Step 3.3.2 -- Plot the Out of Box error by the number of trees 
rf_model2 <- sits_model_export(rf_model)

matplot(rf_model2$err.rate, 
        type = "l", lty = 1, lwd = 2,
        col = my_colors,           
        main = "Out of Box error by the number of trees",
        xlab = "Number of Trees (ntree)", ylab = "Out of Box Error")

legend("topright", 
       legend = names(my_colors), 
       col = my_colors, 
       lty = 1,      
       cex = 1,    
       bty = "n")

ggsave(
  filename = paste0(process_version, "_", tiles_train,"_", no.years, var, "_oob_ntree.png"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# Step 3.4 -- Save the ML model to a R file
saveRDS(rf_model,paste0(rds_path, "model/random_forest/", "RF-model_", length(cube$tile),"-tiles-", tiles_train, "_", no.years,"-period-",cube_dates[1],"_",cube_dates[length(cube_dates)], "_", var, "_", process_version, ".rds"))

print("Model has been trained!")
