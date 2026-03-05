# ============================================================
#  Validation and accuracy measurement
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
model_name       <- "RF-model_4-tiles-012015-012014-013015-013014_2y-period-2023-07-28_2025-07-28_2y-all-classes_2026-02-25_17h58m.rds"
model            <- readRDS(file.path("data/rds/model/random_forest", model_name))
class_dir        <- "data/class"
class_raster_dir <- "data/class/raster"
samples_dir      <- "data/raw/samples/validation_samples"
plots_path       <- "data/plots"
aux_dir          <- "data/raw/auxiliary"
version          <- "rf-2y-012014-all-classes-2y"

dir.create(class_raster_dir, recursive = TRUE, showWarnings = FALSE)

samples_validation_list <- dir(
  samples_dir,
  pattern = paste0(".*_", version, "_.*\\.gpkg$"),
  full.names = TRUE
)

# ============================================================
# 2. Create raster file from classified vector map
# ============================================================

# Step 2.1 -- Function to rasterize
sits_rasterize_segments <- function(file, res, class_raster_dir, style = NULL) {
  
  stopifnot(!is.null(res))
  stopifnot(!is.null(file))
  stopifnot(!is.null(class_raster_dir))
  
  fs::dir_create(class_raster_dir, recurse = TRUE)
  
  file <- fs::path_expand(file)
  class_raster_dir <- fs::path_expand(class_raster_dir)
  
  output_file_base <- fs::path(class_raster_dir) / fs::path_file(file) |>
    fs::path_ext_remove()
  
  output_file <- stringr::str_c(output_file_base, ".tif")
  output_style <- stringr::str_c(output_file_base, ".qml")
  
  if (fs::file_exists(output_file)) {
    return(output_file)
  }
  
  file_sf <- sf::st_read(file, quiet = TRUE)
  
  if (is.null(style)) {
    
    file_sf <- file_sf |>
      dplyr::mutate(
        class_num = .data[["class"]] |>
          as.factor() |>
          as.numeric()
      )
    
    style <- file_sf |>
      tibble::as_tibble() |>
      dplyr::select(dplyr::all_of(c("class", "class_num"))) |>
      dplyr::distinct(.data[["class"]], .data[["class_num"]]) |>
      dplyr::mutate(color = RColorBrewer::brewer.pal(n = dplyr::n(), name = "Set3")) |>
      dplyr::rename(
        index = class_num,
        name = class
      ) |>
      dplyr::arrange(.data[["index"]])
    
  } else {
    
    file_sf <- file_sf |>
      dplyr::rename(name = class) |>
      dplyr::left_join(
        style |> dplyr::select(dplyr::all_of(c("name", "index")))
      ) |>
      dplyr::mutate(
        class_num = .data[["index"]]
      )
  }
  
  file_bbox <- sf::st_bbox(file_sf)
  
  file_gpkg <- fs::file_temp(ext = ".gpkg")
  
  sf::st_write(obj = file_sf, dsn = file_gpkg, quiet = TRUE)
  
  a_srs <- readRDS(
    url("https://github.com/restore-plus/restore-utils/raw/refs/heads/main/inst/extdata/crs/bdc.rds")
  )
  
  gdalUtilities::gdal_rasterize(
    src_datasource = file_gpkg,
    dst_filename = output_file,
    a = "class_num",
    at = TRUE,
    tr = c(res, res),
    te = file_bbox,
    ot = "Int16",
    init = 255,
    a_nodata = 255,
    co = c(
      "COMPRESS=ZSTD",
      "PREDICTOR=2",
      "ZSTD_LEVEL=1",
      "BIGTIFF=YES"
    ),
    a_srs = a_srs
  )
  
  sits:::.colors_qml(
    color_table = style,
    file = output_style
  )
  return(output_file)
}

# Step 2.2 -- Style from ML model
style <- tibble::tibble(
  name = sits_labels(model),
  index = 1:length(sits_labels(model)),
  color = pals::cols25(length(sits_labels(model)))
)

# Step 2.3 -- Rasterize classified vectors
to_raster <- paste0(".*_class_", version, ".*\\.gpkg$")

class_files <- list.files(
  path = class_dir,
  pattern = to_raster,
  full.names = TRUE,
  recursive = TRUE
)

raster_files <- purrr::map(class_files, function(file) {
  file_name <- fs::path_file(file)
  cli::cli_inform("Processing: {file_name}")
  tile_id <- stringr::str_extract(file_name, "\\d{6}")
  period_id <- stringr::str_extract(file_name, "\\d+y")
  tile_period_dir <- file.path(
    class_raster_dir,
    tile_id,
    period_id
  )
  
  fs::dir_create(tile_period_dir, recurse = TRUE)
  sits_rasterize_segments(
    file = file,
    res = 10,
    style = style,
    class_raster_dir = tile_period_dir
  )
})

# ============================================================
# 3. SITS Cube
# ============================================================

# Step 3.1 -- Get labels associated to the trained model data set (Enumerate them in the order they appear according to "sits_labels(model)")
cube_dirs <- list.dirs(class_raster_dir, recursive = TRUE)

cube_dirs <- cube_dirs[
  sapply(cube_dirs, function(x) {
    files <- list.files(x, pattern = "\\.tif$")
    any(grepl(version, files))
  })
]

labels <- c(
  x = sits_labels(model)
)
names(labels) <- 1:length(labels)

# Step 3.2 -- Load the original cube with classified raster file
cube <- sits_cube(
  source = "BDC",
  collection = "SENTINEL-2-16D",
  bands = "class",
  labels = labels,
  data_dir = cube_dirs, # classified raster file cannot be in the same folder as the classified gpkg file
  version = version,
  parse_info = c("satellite", "sensor", "tile", "start_date", "end_date", 
                 "band", "version")
)

# ============================================================
# 4. Full Map Stratified random sampling
# ============================================================

# 4.1 -- Sampling design
sampling_design <- sits_sampling_design(
  cube = cube,
  expected_ua = c(
    "Corpo_Dagua"                           = 0.95,
    "Corte_Raso_Com_Arvores_Remanescentes"  = 0.10,
    "Corte_Raso"                            = 0.70,
    "Corte_Raso_Antigo"                     = 0.85,
    #"Corte_Raso_Com_Fogo"                   = 0.70,
    "Corte_Raso_Com_Vegetacao"              = 0.70,
    "Corte_Raso_Antigo_Com_Vegetacao"       = 0.85,
    "Degradacao"                            = 0.70,
    "Degradacao_Por_Fogo"                   = 0.70,
    "Floresta"                              = 0.95,
    #"Floresta_Transicional"                 = DEFINIR,  
    "Vegetacao_Natural_Nao_Florestal"       = 0.70,
    "Area_Inundavel"                        = 0.70
  ),
  alloc_options = c(120, 100, 75, 50, 30),
  std_err = 0.01,
  rare_class_prop = 0.025
)

# 4.2 -- Show sampling design
sampling_design

# 4.3 -- Generate stratified samples
samples_sf <- sits_stratified_sampling(
  cube = cube,
  sampling_design = sampling_design,
  alloc = "alloc_30",
  overhead = 1.2, # overproportion to avoid border pixels
  progress = TRUE,
  multicores = 12)

# 4.4 -- Total of each class
samples_sf%>% group_by(label) %>% summarise(num = n())

# 4.5 -- Define File Path
samples_sf_file_path <- file.path(samples_dir, paste0("samples-validation-full-map_", version, "_", process_version, ".gpkg"))

# 4.6 -- Save samples_sf object as GPKG file
sf::st_write(samples_sf, samples_sf_file_path, append = FALSE)

# ============================================================
# 5. Accuracy assessment of Full Map classified images
# ============================================================

# Step 5.1 -- Get validation samples points (in geographical coordinates - lat/long)
samples_validation <- st_read(grep("*full-map*", samples_validation_list, value = TRUE)) #full map validation samples

# Step 5.2 -- Calculate accuracy
area_acc_full_map <- sits_accuracy(cube,
                                   validation = samples_validation,
                                   memsize = 180,
                                   multicores = 28) # adapt to your computer CPU core availability

# Step 5.3 -- Print the area estimated accuracy
area_acc_full_map

# Step 5.4 -- Show confusion matrix
area_acc_full_map$error_matrix

# ============================================================
# 6. Plotting Full Map Accuracy
# ============================================================
new_label <-c(
  "Corpo_Dagua"                          = "Water",
  "Corte_Raso_Com_Arvores_Remanescentes" = "Clear Cut Trees", 
  "Corte_Raso"                           = "Clear Cut Bare Soil", 
  "Corte_Raso_Antigo"                    = "Old Clear Cut Bare Soil",
  "Corte_Raso_Com_Fogo"                  = "Fire Suppression", 
  "Corte_Raso_Com_Vegetacao"             = "Clear Cut Vegetation",
  "Corte_Raso_Antigo_Com_Vegetacao"      = "Old Clear Cut Vegetation",
  "Degradacao"                           = "Degradation", 
  "Degradacao_Por_Fogo"                  = "Fire Degradation",
  "Floresta"                             = "Forest",
  "Vegetacao_Natural_Nao_Florestal"      = "Non Forest Natural Vegetation",
  "Floresta_Transicional"                = "Transitional Forest",
  "Area_Inundavel"                       = "Wetland"
)

# To keep the original names
# new_label <- labels

# Step 6.1 -- Create a tibble from error matrix
matriz_conf_full_map <- tibble(as.data.frame(area_acc_full_map$error_matrix))

# Step 6.2 -- Plot error matrix
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
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# Step 6.3 -- Convert accuracies results to a data frame
acuracias_full_map <- data.frame(class = names(area_acc_full_map$accuracy[[1]]),
                                 user_accuracy = round(as.numeric(area_acc_full_map$accuracy[[1]]), 2),
                                 prod_accuracy = round(as.numeric(area_acc_full_map$accuracy[[2]]), 2)
) %>% 
  tidyr::pivot_longer(cols = -class,
                      names_to = "tipo_acuracia",
                      values_to = "acuracia")

# Step 6.4 -- Plot accuracies metrics
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
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# Step 6.5 -- Convert Error-Adjusted Area (ha) results to a data frame
class_areas_full_map <- data.frame(class = names(area_acc_full_map$area_pixels),
                                   mapped_area_ha = round(as.numeric(area_acc_full_map$area_pixels), 2),
                                   error_adj_area_ha = round(as.numeric(area_acc_full_map$error_ajusted_area), 2),
                                   conf_interval_ha = round(as.numeric(area_acc_full_map$conf_interval), 2)
) %>%
  tidyr::pivot_longer(cols = -class,
                      names_to = "tipo_area",
                      values_to = "area")

# Step 6.6 -- Plot Error-Adjusted Area (ha)
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
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# ============================================================
# 7. PRODES Adjusted Map Accuracy
# ============================================================

# 7.1 -- Reclassify classified cube
mask_label <- c("1" = "Natural Vegetation",
                "0" = "Deforestation Mask")

prodes_mask <- sits_cube(source = "BDC",
                         collection = "SENTINEL-2-16D",
                         data_dir = aux_dir,
                         parse_info = c("X1", "X2", "tile", "start_date", "end_date", "band", "version"),
                         bands = "class",
                         version = "v2024",
                         labels = mask_label)

cube_reclass <- purrr::map(cube_dirs_filtered, function(dir_path) {
  
  tile_id <- basename(dirname(dir_path))
  period_id <- basename(dir_path)
  
  cli::cli_inform("Reclassifying tile {tile_id} period {period_id}")
  
  sits_reclassify(
    cube = cube,
    mask = prodes_mask,
    rules = list(
      "Deforestation" =
        cube %in% c(
          "Corte_Raso_Com_Arvores_Remanescentes",
          "Corte_Raso",
          #"Corte_Raso_Com_Fogo",
          "Corte_Raso_Com_Vegetacao"
        ),
      "Other_Classes" =
        cube %in% c(
          "Corpo_Dagua",
          "Corte_Raso_Antigo",
          "Corte_Raso_Antigo_Com_Vegetacao",
          "Floresta",
          "Vegetacao_Natural_Nao_Florestal",
          "Floresta_Transicional",
          "Area_Inundavel",
          "Degradacao",
          "Degradacao_Por_Fogo"
        )
    ),
    multicores = 24,
    memsize = 180,
    version = paste("prodes", version, sep = "-"),
    output_dir = dir_path,
    progress = TRUE
  )
})

# extrair o cube da lista
cube_reclass <- cube_reclass[[1]]

# 7.2 -- Sampling design prodes
sampling_design <- sits_sampling_design(
  cube = cube_reclass,
  expected_ua = c(
    "Deforestation" = 0.70,
    "Other_Classes" = 0.90
  ),
  alloc_options = c(120, 100),
  std_err = 0.01,
  rare_class_prop = 0.05
)

# 7.2 -- Show sampling design
sampling_design

# 7.3 -- Generate stratified samples
samples_sf <- sits_stratified_sampling(
  cube = cube_reclass,
  sampling_design = sampling_design,
  alloc = "alloc_120",
  overhead = 1.2, # overproportion to avoid border pixels
  progress = TRUE,
  multicores = 28)

# 7.4 -- Total of each class
samples_sf%>% group_by(label) %>% summarise(num = n())

# 7.5 -- Define File Path
samples_sf_file_path <- file.path(samples_dir, paste0("samples-validation-prodes-adjusted_", version, "_", process_version, ".gpkg"))

# 7.6 -- Save samples_sf object as GPKG file
sf::st_write(samples_sf, samples_sf_file_path, append = FALSE)

# ============================================================
# 8. Accuracy assessment of PRODES Adjusted Map classified images
# ============================================================

# Step 8.1 -- Get validation samples points (in geographical coordinates - lat/long)
samples_validation <- st_read(grep("*prodes-adjusted*", samples_validation_list, value = TRUE)) #prodes adjusted validation samples

# Step 8.2 -- Calculate accuracy
area_acc_prodes <- sits_accuracy(cube_reclass, 
                                 validation = samples_validation,
                                 memsize = 180,
                                 multicores = 28) # adapt to your computer CPU core availability

# Step 8.3 -- Print the area estimated accuracy
area_acc_prodes

# Step 8.4 -- Show confusion matrix
area_acc_prodes$error_matrix

# ============================================================
# 9. Plotting PRODES Adjusted Map Accuracy
# ============================================================

# Step 9.1 -- Create a tibble from error matrix
matriz_conf_prodes <- tibble(as.data.frame(area_acc_prodes$error_matrix))

# Step 9.2 -- Plot error matrix
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
                   labels = new_label) +
  scale_x_discrete(labels = new_label) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 25,
                                   hjust = 1,
                                   size = 8), 
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5))

ggsave(
  filename = paste(version, "matrix-confusion-prodes.png", sep = "_"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# Step 9.3 -- Convert accuracies results to a data frame
acuracias_prodes <- data.frame(class = names(area_acc_prodes$accuracy[[1]]),
                               user_accuracy = round(as.numeric(area_acc_prodes$accuracy[[1]]), 2),
                               prod_accuracy = round(as.numeric(area_acc_prodes$accuracy[[2]]), 2)
) %>% 
  tidyr::pivot_longer(cols = -class,
                      names_to = "tipo_acuracia",
                      values_to = "acuracia")

# Step 9.4 -- Plot accuracies metrics
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
                   labels = new_label) +
  scale_x_discrete(labels = c("prod_accuracy" = "Prod Acc",
                              "user_accuracy" = "User Acc")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 1, size = 11, margin = margin(t = 10)),
        panel.grid = element_blank())

ggsave(
  filename = paste(version, "accuracies-prodes.png", sep = "_"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# Step 9.5 -- Convert Error-Adjusted Area (ha) results to a data frame
class_areas <- data.frame(class = names(area_acc_prodes$area_pixels),
                          mapped_area_ha = round(as.numeric(area_acc_prodes$area_pixels), 2),
                          error_adj_area_ha = round(as.numeric(area_acc_prodes$error_ajusted_area), 2),
                          conf_interval_ha = round(as.numeric(area_acc_prodes$conf_interval), 2)
) %>%
  tidyr::pivot_longer(cols = -class,
                      names_to = "tipo_area",
                      values_to = "area")

# Step 9.6 -- Plot Error-Adjusted Area (ha)
ggplot(class_areas, aes(x = tipo_area, y = class, fill = area)) + 
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
  filename = paste(version, "areas-prodes.png", sep = "_"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# ============================================================
# 10. PRODES Degradation Adjusted Map Accuracy
# ============================================================

# 10.1 -- Detect tiles and period automatically
tile_version <- stringr::str_extract(version, "\\d{6}")
period_version <- stringr::str_extract(version, "\\d+y")

cube_dirs_filtered <- cube_dirs[
  grepl(tile_version, cube_dirs) &
    grepl(period_version, cube_dirs)
]

cube_reclass <- purrr::map(cube_dirs_filtered, function(dir_path) {
  tile_id <- basename(dirname(dir_path))
  period_id <- basename(dir_path)
  cli::cli_inform("Reclassifying tile {tile_id} period {period_id}")
  sits_reclassify(
    cube = cube,
    mask = prodes_mask,
    rules = list(
      "Deforestation" =
        cube %in% c(
          "Corte_Raso_Com_Arvores_Remanescentes",
          "Corte_Raso",
          "Corte_Raso_Com_Vegetacao"
        ),
      "Degradation" =
        cube %in% c(
          "Degradacao",
          "Degradacao_Por_Fogo"
        ),
      "Other_Classes" =
        cube %in% c(
          "Corpo_Dagua",
          "Corte_Raso_Antigo",
          "Corte_Raso_Antigo_Com_Vegetacao",
          "Floresta",
          "Vegetacao_Natural_Nao_Florestal",
          "Floresta_Transicional",
          "Area_Inundavel"
        )
    ),
    multicores = 24,
    memsize = 180,
    version = paste("prodes-degradation", version, sep = "-"),
    output_dir = dir_path,
    progress = TRUE
  )
})

# 10.2 -- Extract the cube REQUIRED
cube_reclass <- cube_reclass[[1]]

# 10.2 -- Sampling design degradation
sampling_design <- sits_sampling_design(
  cube = cube_reclass,
  expected_ua = c(
    "Deforestation" = 0.70,
    "Degradation" = 0.60, 
    "Other_Classes" = 0.95
  ),
  alloc_options = c(120, 100),
  std_err = 0.01,
  rare_class_prop = 0.05
)

# 10.2 -- Show sampling design
sampling_design

# 10.3 -- Generate stratified samples
samples_sf <- sits_stratified_sampling(
  cube = cube_reclass,
  sampling_design = sampling_design,
  alloc = "alloc_100",
  overhead = 1.2, # overproportion to avoid border pixels
  progress = TRUE,
  multicores = 28)

# 10.4 -- Total of each class
samples_sf%>% group_by(label) %>% summarise(num = n())

# 10.5 -- Define File Path
samples_sf_file_path <- file.path(samples_dir, paste0("samples-validation-desmat-degrad_", version, "_", process_version, ".gpkg"))

# 10.6 -- Save samples_sf object as GPKG file
sf::st_write(samples_sf, samples_sf_file_path, delete_dsn = TRUE, append = FALSE)

# ============================================================
# 11. Accuracy assessment of PRODES Degradation Adjusted Map classified images
# ============================================================

# Step 11.1 -- Get validation samples points (in geographical coordinates - lat/long)
samples_validation <- st_read(grep("*prodes-degrad*", samples_validation_list, value = TRUE)) #prodes adjusted validation samples with degradation classes

# Step 11.2 -- Calculate accuracy
area_acc_prodes <- sits_accuracy(cube_reclass, 
                                 validation = samples_validation,
                                 memsize = 180,
                                 multicores = 28) # adapt to your computer CPU core availability

# Step 11.3 -- Print the area estimated accuracy
area_acc_prodes

# Step 11.4 -- Show confusion matrix
area_acc_prodes$error_matrix

# ============================================================
# 12. Plotting PRODES Degradation Adjusted Map Accuracy
# ============================================================

# Step 12.1 -- Create a tibble from error matrix
matriz_conf_prodes <- tibble(as.data.frame(area_acc_prodes$error_matrix))

# Step 12.2 -- Plot error matrix
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
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# Step 12.3 -- Convert accuracies results to a data frame
acuracias_prodes <- data.frame(class = names(area_acc_prodes$accuracy[[1]]),
                               user_accuracy = round(as.numeric(area_acc_prodes$accuracy[[1]]), 2),
                               prod_accuracy = round(as.numeric(area_acc_prodes$accuracy[[2]]), 2)
) %>% 
  tidyr::pivot_longer(cols = -class,
                      names_to = "tipo_acuracia",
                      values_to = "acuracia")

# Step 12.4 -- Plot accuracies metrics
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
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# Step 12.5 -- Convert Error-Adjusted Area (ha) results to a data frame 
class_areas <- data.frame(class = names(area_acc_prodes$area_pixels),
                          mapped_area_ha = round(as.numeric(area_acc_prodes$area_pixels), 2),
                          error_adj_area_ha = round(as.numeric(area_acc_prodes$error_ajusted_area), 2),
                          conf_interval_ha = round(as.numeric(area_acc_prodes$conf_interval), 2)
) %>%
  tidyr::pivot_longer(cols = -class,
                      names_to = "tipo_area",
                      values_to = "area")

# Step 12.6 -- Plot Error-Adjusted Area (ha)
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
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)