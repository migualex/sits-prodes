# ============================================================
# Train Machine Learning Model
# ============================================================

## I. Load Required Libraries
library(sits)
library(ggplot2)

## II. Define the date and time for the start of processing
date_process <- format(Sys.Date(), "/%Y_%m_%d_")
time_process <- format(Sys.time(), "%Hh%Mm_", tz = "America/Sao_Paulo")
process_version <- paste0(date_process, time_process)

## III. Define the paths for files and folders needed in the processing
rds_path      <- "data/rds"
plots_path    <- "data/plots"

# IV. Identifier to distinguish this model run from previous versions | Keep it the same as the samples' identifier
var <- "with-df-mask"

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


# ============================================================
# 2. Training and saving model
# ============================================================

# Step 2.1 -- Read RDS file that contains clean and balanced training samples
clean_samples_balanced <- readRDS(
  paste0(rds_path, "/time_series/", process_version, tiles_train, var, "_clean_samples_balanced", ".rds")
)

# Step 2.2 -- Set a seed of random number generator (RNG) for reproducibility
set.seed(88)

# Step 2.3 -- Train the model
rf_model <- sits_train(
   samples   = clean_samples_balanced,
   ml_method = sits_rfor(num_trees = 100)
 )

# Step 2.3.1 -- Plot the most important variables of the model
plot(rf_model)

# Step 2.3.2 -- Plot the Out of Box error by the number of trees 
rf_model2 <- sits_model_export(rf_model)

matplot(rf_model2$err.rate, 
        type = "l", lty = 1, lwd = 2,
        col = my_colours,           
        main = "Out of Box error by the number of trees",
        xlab = "Number of Trees (ntree)", ylab = "Out of Box Error")

legend("topright", 
       legend = names(my_colours), 
       col = my_colours, 
       lty = 1,      
       cex = 1,    
       bty = "n")

ggsave(
  filename = paste0(process_version, tiles_train,"_", no.years, var, "_oob_ntree.png"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# Step 2.4 -- Save the ML model to a R file
saveRDS(rf_model,paste0(rds_path, "/model/random_forest/", process_version, "RF-", no.years, "-", tiles_train, "-",  var,".rds"))

print("Model has been trained!")