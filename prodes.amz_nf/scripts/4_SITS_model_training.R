# ============================================================
# Train Machine Learning Model
# ============================================================

# ============================================================
# 1. Libraries, paths and some initial parameters
# ============================================================

# Step 1.1 -- Load Required Libraries
library(sits)
library(ggplot2)

# Step 1.2 -- Define the date and time for the start of processing
date_process <- format(Sys.Date(), "%Y-%m-%d_")
time_process <- format(Sys.time(), "%Hh%Mm", tz = "America/Sao_Paulo")
process_version <- paste0(date_process, time_process)

# Step 1.3 -- Define the paths for files and folders needed in the processing
time_series_name  <- "samples_2-tiles-014002-015002_2y-period-2023-07-28_2025-07-28_nf_2026-03-03_15h08m.rds" #add the time series name
rds_path          <- "data/rds/"
plots_path        <- "data/plots/"
time_series_path  <- file.path("data/rds/time_series/", time_series_name)

# Step 1.4 -- Identifier to distinguish this model run from previous versions | Keep it the same as the samples' identifier
var <- "nf-samples-crude"

# Step 1.5 -- Define a list with preference colours for each class
my_colors <- c(
  "OOB"                                                                            = "black",
  "Hidrografia_Rio"                                                                = "#1f78b4",
  "Hidrografia_Lago"                                                               = "#2980B9",
  "Vegetacao_Natural_Nao_Florestal_Herbacea_Umida"                                 = "#A0B9C8",
  "Vegetacao_Natural_Nao_Florestal_Herbacea_Seca_Mais_Biomassa"                    = "#f6cc41",
  "Vegetacao_Natural_Nao_Florestal_Herbacea_Seca_Menos_Biomassa"                   = "#DBEBD8",
  "Vegetacao_Natural_Nao_Florestal_Vereda"                                         = "#3ABABA",
  "Vegetacao_Natural_Nao_Florestal_Mata"                                           = "#1E8449",
  "Fogo_Recente_Em_Vegetacao_Natural_Nao_Florestal"                                = "#CD6155",
  "Fogo_Antigo_Em_Vegetacao_Natural_Nao_Florestal"                                 = "#E6B0AA",
  "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Agricultura"                   = "#F0B27A",
  "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Reservatorio"                  = "#6D98B8",
  "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Solo_Exposto"                  = "#F39C12",
  "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Agricultura_Antigo"            = "#B08B57",
  "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Reservatorio_Antigo"           = "#4A677D",
  "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Solo_Exposto_Antigo"           = "#A0522D"
)

# Step 1.6 -- Define time range
start_date    <- "2023-08-01"
end_date      <- "2025-07-31"

# ============================================================
# 2. Define and Load Data Cubes
# ============================================================

# Step 2.1 -- Create a training cube from a collection
cube <- sits_cube(
  source      = "BDC",
  collection  = "SENTINEL-2-16D",
  bands       = c('B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12', 'NDVI', 'NBR', 'EVI', 'CLOUD'),
  tiles       = c('014002', '015002'),
  start_date  = start_date,
  end_date    = end_date,
  progress    = TRUE
)

# Step 2.2 -- Calculate the number of years in the training cube
cube_dates <- sits_timeline(cube)
no.years <- paste0(floor(lubridate::year(end_date) - lubridate::year(start_date)), "y")

# Step 2.3 -- Concatenates all the names of the training tiles into a single string separated by '-'
tiles_train <- paste(cube$tile, collapse = "-")

# 2.4 Create output directory per tile and period
tiles_id <- paste(sort(unique(tiles_train)), collapse = "_")

tile_period_dir <- file.path(plots_path, tiles_id, no.years)

dir.create(tile_period_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 3. Cross-validation of training data
# ============================================================

# Step 3.1 -- Reading training samples
train_samples <- readRDS(time_series_path)

# Step 3.2 -- Using k-fold validation
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

# Step 3.2.1 -- Plot the confusion matrix
plot(rfor_validate, type = "confusion_matrix")

# Step 3.2.2 -- Plot the metrics by class
plot(rfor_validate, type = "metrics")

# ============================================================
# 4. Training and saving model
# ============================================================

# Step 4.1 -- Set a seed of random number generator (RNG) for reproducibility
set.seed(88)

# Step 4.2 -- Train the model
rf_model <- sits_train(
  samples   = train_samples,
  ml_method = sits_rfor(num_trees = 100)
)

# Step 4.2.1 -- Plot the most important variables of the model
plot(rf_model)

# Step 4.2.2 -- Save the plot
ggsave(
  filename = paste0(process_version, "_", tiles_train,"_", no.years, var, "_minimal_tree_depth_mde.png"),
  path = tile_period_dir,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# Step 4.3 --  Exports the model as an object for further exploration
rf_model2 <- sits_model_export(rf_model)

# Step 4.3.1 -- Save the plot
png(
  filename = file.path(
    tile_period_dir,
    paste0(process_version, "_", tiles_train, "_", no.years, var, "_oob_ntree_mde.png")
  ),
  width = 3529,
  height = 1578,
  res = 350
)

# Step 4.3.2 -- Plot the Out of Box error by the number of trees
matplot(rf_model2$err.rate, 
        type = "l", lty = 1, lwd = 2,
        col = my_colors,           
        main = "Out of Box error by the number of trees",
        xlab = "Number of Trees (ntree)", 
        ylab = "Out of Box Error")

# Step 4.3.3 -- Adding legend to plot
legend("topright", 
       legend = names(my_colors), 
       col = my_colors, 
       lty = 1,      
       cex = 1,    
       bty = "n")
dev.off()

# Step 4.4 -- Save the ML model to a R file
saveRDS(rf_model,paste0(rds_path, "model/random_forest/", "RF-model_", length(cube$tile),"-tiles-", tiles_train, "_", no.years,"-period-",cube_dates[1],"_",cube_dates[length(cube_dates)], "_", var, "_", process_version, ".rds"))

print("Model has been trained!")