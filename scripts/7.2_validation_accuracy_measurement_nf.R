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
library(fs)
library(stringr)
library(purrr)

# Step 1.2 -- Define the date and time for the start of processing
date_process <- format(Sys.Date(), "%Y-%m-%d_")
time_process <- format(Sys.time(), "%Hh%Mm", tz = "America/Sao_Paulo")
process_version <- paste0(date_process, time_process)

# Step 1.3 -- Define the paths for files and folders needed in the processing
model_name       <- "RF-model_2-tiles-014002-015002_2y-period-2023-07-28_2025-07-28_nf-samples-crude_2026-03-05_08h56m.rds"
model            <- readRDS(file.path("data/rds/model/random_forest", model_name))
class_dir        <- "data/class"
samples_dir      <- "data/raw/samples/validation_samples"
plots_path       <- "data/plots"
aux_dir          <- "data/raw/auxiliary"
version          <- "prodes-degradation-rf-2y-novos-segmentos"
plots_dir        <- file.path(plots_path, version)
class_raster_dir <- file.path(class_dir, str_split_i(version, pattern = "-", 3), "raster")

# Step 1.4 -- Create the directory for storing class rasters, including any necessary parent directories. Suppress warnings if the directory already exists.
dir.create(class_raster_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# Step 1.5 -- Get the list of validation sample files matching the version pattern in the samples directory
samples_validation_list <- dir(
  samples_dir,
  pattern = paste0(".*", str_split_i(version, pattern = "-", 3), ".*\\.gpkg$"),
  full.names = TRUE
)

# ============================================================
# 2. SITS Cube
# ============================================================

# Step 2.1 -- Get labels associated to the trained model data set (Enumerate them in the order they appear according to "sits_labels(model)")
cube_dirs <- list.dirs(class_raster_dir, recursive = TRUE)

# Step 2.2 -- # Filters the list of directories, keeping only those that contain .tif files matching the specified version.
cube_dirs <- cube_dirs[
    sapply(cube_dirs, function(x) {
      files <- list.files(x, pattern = "\.tif$")
      any(grepl(version, files))
    })
]

# Step 2.3 -- Store the labels from the trained model using sits_labels and assign them to a named vector
labels <- c(
  x = sits_labels(model)
)
names(labels) <- 1:length(labels)

# Step 2.4 -- Load the original cube with classified raster file
cube <- sits_cube(
  source = "BDC",
  collection = "SENTINEL-2-16D",
  bands = "class",
  labels = labels,
  data_dir = "data/class/temp/all_tiles",
  version = version,
  parse_info = c("satellite", "sensor", "tile", "start_date", "end_date", 
                 "band", "version"))

# ============================================================
# 3. Accuracy assessment of Full Map classified images
# ============================================================

# Step 3.1 -- Get validation samples points (in geographical coordinates - lat/long)
samples_validation <- st_read(grep("*full-map*", samples_validation_list, value = TRUE)[1]) #full map validation samples

# Step 3.2 -- Calculate accuracy
area_acc_full_map <- sits_accuracy(cube,
                                   validation = samples_validation,
                                   memsize = 180,
                                   multicores = 28) # adapt to your computer CPU core availability

# Step 3.3 -- Print the area estimated accuracy
area_acc_full_map

# Step 3.4 -- Show confusion matrix
area_acc_full_map$error_matrix

# ============================================================
# 4. Plotting Full Map Accuracy
# ============================================================

# Step 4.1 -- Create a tibble from error matrix
matriz_conf_full_map <- tibble(as.data.frame(area_acc_full_map$error_matrix))

# Step 4.2 -- Plot error matrix
ggplot(matriz_conf_full_map, aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 4) +
  scale_fill_distiller(palette = "Blues",
                       direction = 1,
                       name = "Cases") +
  labs( x = "Reference",
        y = "Predicted", 
        title = "Confusion Matrix") +
  scale_y_discrete(limits = rev,
                   labels = new_label) +
  scale_x_discrete(labels = new_label) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 25,
                                   hjust = 1,
                                   size = 8), 
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5))

ggsave(
  filename = paste(version, "matrix-confusion-full-map.png", sep = "_"),
  path = plots_dir,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# Step 4.3 -- Convert accuracies results to a data frame
acuracias_full_map <- data.frame(class = names(area_acc_full_map$accuracy[[1]]),
                                 user_accuracy = round(as.numeric(area_acc_full_map$accuracy[[1]]), 2),
                                 prod_accuracy = round(as.numeric(area_acc_full_map$accuracy[[2]]), 2)
) %>% 
  tidyr::pivot_longer(cols = -class,
                      names_to = "tipo_acuracia",
                      values_to = "acuracia")

# Step 4.4 -- Plot accuracies metrics
ggplot(acuracias_full_map, aes(x = tipo_acuracia, y = class, fill = acuracia)) + 
  geom_tile(color = "white") +
  geom_text(aes(label = acuracia), size = 4) +
  scale_fill_distiller(palette = "RdYlGn",
                       direction = 1,
                       name = "Value") +
  labs(y = "Class",
       x = "Accuracy", 
       title = "Accuracies",
       caption = paste0("Global Accuracy: ", round(area_acc_full_map$accuracy[[3]], 2))) +
  scale_y_discrete(limits = rev,
                   labels = new_label) +
  scale_x_discrete(labels = c("prod_accuracy" = "Prod Acc",
                              "user_accuracy" = "User Acc")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 1, size = 11, margin = margin(t = 10)),
        panel.grid = element_blank())

ggsave(
  filename = paste(version, "accuracies-full-map.png", sep = "_"),
  path = plots_dir,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# Step 4.5 -- Convert Error-Adjusted Area (ha) results to a data frame
class_areas_full_map <- data.frame(class = names(area_acc_full_map$area_pixels),
                                   mapped_area_ha = round(as.numeric(area_acc_full_map$area_pixels), 2),
                                   error_adj_area_ha = round(as.numeric(area_acc_full_map$error_ajusted_area), 2),
                                   conf_interval_ha = round(as.numeric(area_acc_full_map$conf_interval), 2)
) %>%
  tidyr::pivot_longer(cols = -class,
                      names_to = "tipo_area",
                      values_to = "area")

# Step 4.6 -- Plot Error-Adjusted Area (ha)
ggplot(class_areas_full_map, aes(x = tipo_area, y = class, fill = area)) + 
  geom_tile(color = NA, fill = "white") +
  geom_text(aes(label = area), size = 4) +
  geom_vline(xintercept = c(1.5, 2.5), color = "grey30") +
  labs(y = "Class",
       x = "Metrics", 
       title = "Area Metrics") +
  scale_y_discrete(limits = rev,
                   labels = new_label) +
  scale_x_discrete(limits = rev,
                   labels = c("mapped_area_ha" = "Mapped Area (ha)",
                              "error_adj_area_ha" = "Error-Adjusted Area (ha)",
                              "conf_interval_ha" = "Conf Interval (ha)")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        panel.grid = element_blank())

ggsave(
  filename = paste(version, "areas-full-map.png", sep = "_"),
  path = plots_dir,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# ============================================================
# 5. PRODES Degradation Adjusted Map Accuracy
# ============================================================
# 5.1 -- Reclassify classified cube
mask_label <- c("1" = "Natural Vegetation",
                "0" = "Deforestation Mask")

prodes_mask <- sits_cube(source = "BDC",
                         collection = "SENTINEL-2-16D",
                         data_dir = aux_dir,
                         parse_info = c("X1", "X2", "tile", "start_date", "end_date", "band", "version"),
                         bands = "class",
                         version = "v2024",
                         labels = mask_label)

# 5.2 -- Detect tiles and period automatically
tile_version <- stringr::str_extract(version, "\\d{6}")
period_version <- stringr::str_extract(version, "\\d+y")

cube_dirs_filtered <- cube_dirs[
  grepl(tile_version, cube_dirs)
]

# File name 
class_cube <- sits_cube(
  source = "BDC",
  collection = "SENTINEL-2-16D",
  bands = "class",
  labels = c("15" = "Deforestation", "16" = "Water", 
             "17" = "Grassland", "18" =  "Forest"),
  tiles =  "000000",
  start_date = "2023-08-12",
  end_date = "2025-07-28",
  version = "prodes-degradation-rf-2y-novos-segmentos",
  data_dir = "~/grupos/biomasbr/amazonia/sits-prodes/prodes.amz_nf/data/class/temp/all_tiles",
  parse_info = c("satellite", "sensor", "tile", "start_date", "end_date", 
                 "band", "version"))

# ============================================================
# 6. Accuracy assessment of PRODES Degradation Adjusted Map classified images
# ============================================================

# Step 6.1 -- Get validation samples points (in geographical coordinates - lat/long)
samples_validation <- st_read(grep("*7-tiles-desmat-degrad*", samples_validation_list, value = TRUE))

# Step 6.2 -- Calculate accuracy
area_acc_prodes <- sits_accuracy(class_cube, 
                                 validation = samples_validation,
                                 memsize = 180,
                                 multicores = 28) # adapt to your computer CPU core availability

# Step 6.3 -- Print the area estimated accuracy
area_acc_prodes

# Step 6.4 -- Show confusion matrix
area_acc_prodes$error_matrix

# ============================================================
# 7. Plotting PRODES Degradation Adjusted Map Accuracy
# ============================================================

# Step 7.1 -- Create a tibble from error matrix
matriz_conf_prodes <- tibble(as.data.frame(area_acc_prodes$error_matrix))

# Step 7.2 -- Plot error matrix
ggplot(matriz_conf_prodes, aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 4) +
  scale_fill_distiller(palette = "Blues",
                       direction = 1,
                       name = "Cases") +
  labs( x = "Reference",
        y = "Predicted", 
        title = "Confusion Matrix") +
  scale_y_discrete(limits = rev,
                   labels = c("Other_Classes" = "Other Classes")) +
  scale_x_discrete(labels = c("Other_Classes" = "Other Classes")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 25,
                                   hjust = 1,
                                   size = 8), 
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5))

ggsave(
  filename = paste(version, "matrix-confusion-prodes-degrad.png", sep = "_"),
  path = plots_dir,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# Step 7.3 -- Convert accuracies results to a data frame
acuracias_prodes <- data.frame(class = names(area_acc_prodes$accuracy[[1]]),
                               user_accuracy = round(as.numeric(area_acc_prodes$accuracy[[1]]), 2),
                               prod_accuracy = round(as.numeric(area_acc_prodes$accuracy[[2]]), 2)
) %>% 
  tidyr::pivot_longer(cols = -class,
                      names_to = "tipo_acuracia",
                      values_to = "acuracia")

# Step 7.4 -- Plot accuracies metrics
ggplot(acuracias_prodes, aes(x = tipo_acuracia, y = class, fill = acuracia)) + 
  geom_tile(color = "white") +
  geom_text(aes(label = acuracia), size = 4) +
  scale_fill_distiller(palette = "RdYlGn",
                       direction = 1,
                       name = "Value") +
  labs(y = "Class",
       x = "Accuracy", 
       title = "Accuracies",
       caption = paste0("Global Accuracy: ", round(area_acc_prodes$accuracy[[3]], 2))) +
  scale_y_discrete(limits = rev,
                   labels = c("Other_Classes" = "Other Classes")) +
  scale_x_discrete(labels = c("prod_accuracy" = "Prod Acc",
                              "user_accuracy" = "User Acc")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 1, size = 11, margin = margin(t = 10)),
        panel.grid = element_blank())

ggsave(
  filename = paste(version, "accuracies-prodes-degrad.png", sep = "_"),
  path = plots_dir,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# Step 7.5 -- Convert Error-Adjusted Area (ha) results to a data frame 
class_areas <- data.frame(class = names(area_acc_prodes$area_pixels),
                          mapped_area_ha = round(as.numeric(area_acc_prodes$area_pixels), 2),
                          error_adj_area_ha = round(as.numeric(area_acc_prodes$error_ajusted_area), 2),
                          conf_interval_ha = round(as.numeric(area_acc_prodes$conf_interval), 2)
) %>%
  tidyr::pivot_longer(cols = -class,
                      names_to = "tipo_area",
                      values_to = "area")

# Step 7.6 -- Plot Error-Adjusted Area (ha)
ggplot(class_areas, aes(x = tipo_area, y = class, fill = area)) + 
  geom_tile(color = NA, fill = "white") +
  geom_text(aes(label = area), size = 4) +
  geom_vline(xintercept = c(1.5, 2.5), color = "grey30") +
  labs(y = "Class",
       x = "Metrics", 
       title = "Area Metrics") +
  scale_y_discrete(limits = rev,
                   labels = c("Other_Classes" = "Other Classes")) +
  scale_x_discrete(limits = rev,
                   labels = c("mapped_area_ha" = "Mapped Area (ha)",
                              "error_adj_area_ha" = "Error-Adjusted Area (ha)",
                              "conf_interval_ha" = "Conf Interval (ha)")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        panel.grid = element_blank())

ggsave(
  filename = paste(version, "areas-prodes-degrad.png", sep = "_"),
  path = plots_dir,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)