# ============================================================
#  Validation and accuracy measurement in Non-Forest Areas
# ============================================================

# ============================================================
# 1. Libraries, paths and some initial parameters
# ============================================================

# Step 1.1 -- Load Required Libraries
library(tibble)
library(sits)
library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(stringr)

# Step 1.2 -- Define the paths for files and folders needed in the processing
model_name       <- "model-name.rds"
model            <- readRDS(file.path("data/rds/model/random_forest", model_name))
class_dir        <- "data/class"
samples_dir      <- "data/raw/samples/validation_samples/012014"
plots_dir       <- "data/plots"
mask_dir          <- "data/raw/auxiliary/masks"
version          <- "version"

# Step 1.3 -- Get the list of validation sample files matching the version pattern in the samples directory
samples_validation_list <- dir(
  samples_dir,
  pattern = paste0(".*", str_split_i(version, pattern = "-", 3), ".*\\.gpkg$"),
  full.names = TRUE
)

# Step 1.4 -- Define the list of tiles and period
tiles = c('012014')
start_date = "2023-08-13"
end_date = "2025-07-28"

# Step 1.5 -- Define plotting function
plot_accuracy <- function(acc, version, tile, plots_dir, prefix) {
  
  today     <- format(Sys.Date(), "%Y-%m-%d")
  save_plot <- function(filename, width = 3529, height = 1578) {
    ggsave(
      filename = paste0(filename, "-", prefix, "_", tile, "_", version, "_", today, ".png"),
      path     = plots_dir,
      scale    = 1,
      width    = width,
      height   = height,
      units    = "px",
      dpi      = 350
    )
  }
  
  label_other <- c("Other_Classes" = "Other Classes")
  
  # ── 1. Confusion Matrix ───────────────────────────────────────────────────────
  matriz_conf <- tibble(as.data.frame(acc$error_matrix))
  
  p1 <- ggplot(matriz_conf, aes(x = Var2, y = Var1, fill = Freq)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Freq), size = 4) +
    scale_fill_distiller(palette = "Blues", direction = 1, name = "Cases") +
    labs(x = "Reference", y = "Predicted", title = "Confusion Matrix") +
    scale_y_discrete(limits = rev, labels = label_other) +
    scale_x_discrete(labels = label_other) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      plot.title  = element_text(hjust = 0.5)
    )
  
  print(p1)
  save_plot("confusion-matrix")
  
  # ── 2. Accuracies ─────────────────────────────────────────────────────────────
  acuracias <- data.frame(
    class         = names(acc$accuracy[[1]]),
    user_accuracy = round(as.numeric(acc$accuracy[[1]]), 2),
    prod_accuracy = round(as.numeric(acc$accuracy[[2]]), 2)
  ) %>%
    dplyr::mutate(
      f1_score = round(2 * (user_accuracy * prod_accuracy) / (user_accuracy + prod_accuracy), 2)
    ) %>%
    tidyr::pivot_longer(cols = -class, names_to = "tipo_acuracia", values_to = "acuracia")
  
  p2 <- ggplot(acuracias, aes(x = tipo_acuracia, y = class, fill = acuracia)) +
    geom_tile(color = "white") +
    geom_text(aes(label = acuracia), size = 4) +
    scale_fill_distiller(palette   = "RdYlGn",
                         direction = 1,
                         name      = "Value",
                         limits    = c(0, 1)) +
    labs(
      y       = "Class",
      x       = "Accuracy",
      title   = "Accuracies",
      caption = paste0("Global Accuracy: ", round(acc$accuracy[[3]], 2))
    ) +
    scale_y_discrete(limits = rev, labels = label_other) +
    scale_x_discrete(labels = c("prod_accuracy" = "Prod Acc",
                                "user_accuracy" = "User Acc",
                                "f1_score"      = "F1 Score")) +
    theme_minimal() +
    theme(
      plot.title   = element_text(hjust = 0.5),
      plot.caption = element_text(hjust = 1, size = 11, margin = margin(t = 10)),
      panel.grid   = element_blank()
    )
  
  print(p2)
  save_plot("metrics")
  
  # ── 3. Area Metrics ───────────────────────────────────────────────────────────
  class_areas <- data.frame(
    class             = names(acc$area_pixels),
    mapped_area_ha    = round(as.numeric(acc$area_pixels), 2),
    error_adj_area_ha = round(as.numeric(acc$error_ajusted_area), 2),
    conf_interval_ha  = round(as.numeric(acc$conf_interval), 2)
  ) %>%
    tidyr::pivot_longer(cols = -class, names_to = "tipo_area", values_to = "area")
  
  p3 <- ggplot(class_areas, aes(x = tipo_area, y = class, fill = area)) +
    geom_tile(color = NA, fill = "white") +
    geom_text(aes(label = area), size = 4) +
    geom_vline(xintercept = c(1.5, 2.5), color = "grey30") +
    labs(y = "Class", x = "Metrics", title = "Area Metrics") +
    scale_y_discrete(limits = rev, labels = label_other) +
    scale_x_discrete(
      limits = rev,
      labels = c(
        "mapped_area_ha"    = "Mapped Area (ha)",
        "error_adj_area_ha" = "Error-Adjusted Area (ha)",
        "conf_interval_ha"  = "Conf Interval (ha)"
      )
    ) +
    theme_minimal() +
    theme(
      plot.title      = element_text(hjust = 0.5),
      legend.position = "none",
      panel.grid      = element_blank()
    )
  
  print(p3)
  save_plot("areas")
  
  message("Plots salvos em: ", plots_dir)
  
  invisible(list(confusion_matrix = p1, accuracies = p2, area_metrics = p3))
}

# ============================================================
# 2. Accuracy assessment of Full Map classified images
# ============================================================

# Step 2.1 -- Retrieve local cube of Full Map classified
cube <- sits_cube(
  source = "BDC",
  collection = "SENTINEL-2-16D",
  bands = "class",
  labels = c("1"  = "Wetland", # List the classes according to the number sequence in which they appear in your raster
             "2"  = "Water",
             "3"  = "Forest",
             "4"  = "Transition_Forest",
             "5"  = "Non_Forest_Natural_Vegetation",
             "6"  = "Degradation",
             "7"  = "Degradation_Fire",
             "8"  = "Clear_Cut_Bare_Soil",
             "9"  = "Clear_Cut_Vegetation",
             "10" = "Clear_Cut_Burned_Area",
             "11" = "Clear_Cut_Trees",
             "12" = "Previous_Clear_Cut_Bare_Soil",
             "13" = "Previous_Clear_Cut_Vegetation"
             ),
  tiles =  tiles,
  start_date = start_date,
  end_date = end_date,
  version = "raster-version",
  data_dir = "data/class/012014/original_class/012014/raster",
  parse_info = c("satellite", "sensor", "tile", "start_date", "end_date", 
                 "band", "version"))

# Step 2.2 -- Get validation samples points (in geographical coordinates - lat/long)
samples_validation <- st_read(grep("*samples-version*", samples_validation_list, value = TRUE)[1]) #full map validation samples

# Step 2.3 -- Calculate accuracy
full_map_acc <- sits_accuracy(cube,
                              validation = samples_validation,
                              memsize = 180,
                              multicores = 28) # adapt to your computer CPU core availability

# Step 2.4 -- Print the area estimated accuracy
full_map_acc

# Step 2.5 -- Show confusion matrix
full_map_acc$error_matrix

# Step 2.6 -- Plotting Full Map Accuracy
plot_accuracy(
  acc       = full_map_acc,
  version   = "model-version",
  tile      = tiles,
  plots_dir = plots_dir,
  prefix    = "full-map-acc"
)

# ============================================================
# 3. Accuracy assessment of PRODES Adjusted Map classified
# ============================================================

# Step 3.1 -- Retrieve local cube of PRODES adjusted map classified
class_cube <- sits_cube(
  source = "BDC",
  collection = "SENTINEL-2-16D",
  bands = "class",
  labels = c("14" = "Deforestation",
             "15" = "Degradation",
             "16" = "Other_Classes",
             ),
  tiles =  tiles,
  start_date = start_date,
  end_date = end_date,
  version = "raster-version",
  data_dir = "data/class/012014/original_class/012014/raster",
  parse_info = c("satellite", "sensor", "tile", "start_date", "end_date", 
                 "band", "version"))

# Step 3.2 -- Get validation samples points (in geographical coordinates - lat/long)
samples_validation <- st_read(grep("*samples-version*", samples_validation_list, value = TRUE))

# Step 3.3 -- Calculate accuracy
prodes_acc <- sits_accuracy(class_cube, 
                            validation = samples_validation,
                            memsize = 180,
                            multicores = 28) # adapt to your computer CPU core availability

# Step 3.4 -- Print the area estimated accuracy
prodes_acc

# Step 3.5 -- Show confusion matrix
prodes_acc$error_matrix

# Step 3.6 -- Plotting PRODES Adjusted Map Accuracy
plot_accuracy(
  acc       = prodes_acc,
  version   = "model-version",
  tile      = tiles,
  plots_dir = plots_dir,
  prefix    = "prodes-acc"
)
