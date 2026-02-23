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

# Step 1.2 -- Define the date and time for the start of processing
date_process <- format(Sys.Date(), "%Y-%m-%d_")
time_process <- format(Sys.time(), "%Hh%Mm", tz = "America/Sao_Paulo")
process_version <- paste0(date_process, time_process)

# Step 1.3 -- Define the paths for files and folders needed in the processing
class_dir <- "data/class"
class_raster_output <- "data/class/raster" # classified raster file cannot be in the same folder as the classified gpkg file
samples_dir <- "data/raw/samples/validation_samples/"
plots_path    <- "data/plots/"
aux_dir <- "data/raw/auxiliary"
model <- readRDS("data/rds/model/random_forest/RF-model_2-tiles-014002-015002_0y-period-2024-07-27_2025-07-28_nf_2026-02-20_14h43m.rds")
version <- "rf-0y-014002-nf-samples-crude"


# ============================================================
# 2. Create raster file from classified vector map
# ============================================================

# Step 2.1 -- Function to rasterize
sits_rasterize_segments <- function(file, res, class_raster_dir, style = NULL) {
  stopifnot(!is.null(res))
  stopifnot(!is.null(file))
  stopifnot(!is.null(class_raster_dir))
  
  # create output dir
  fs::dir_create(class_raster_dir, recurse = TRUE)
  
  # expand paths
  file <- fs::path_expand(file)
  class_raster_dir <- fs::path_expand(class_raster_dir)
  
  # define output files
  output_file_base <- fs::path(class_raster_dir) / fs::path_file(file) |>
    fs::path_ext_remove()
  
  output_file <- stringr::str_c(output_file_base, ".tif")
  output_style <- stringr::str_c(output_file_base, ".qml")
  
  if (fs::file_exists(output_file)) {{
    return(output_file)
  }}
  
  file_sf <- sf::st_read(file, quiet = TRUE)
  
  if (is.null(style)) {
    file_sf <- file_sf |>
      dplyr::mutate(class_num = .data[["class"]] |>
                      as.factor() |>
                      as.numeric())
    
    style <- file_sf |>
      tibble::as_tibble() |>
      dplyr::select(dplyr::all_of(c("class", "class_num"))) |>
      dplyr::distinct(.data[["class"]], .data[["class_num"]]) |>
      dplyr::mutate(color = RColorBrewer::brewer.pal(n = dplyr::n(), name = "Set3")) |>
      dplyr::rename("index" = "class_num",
                    "color" = "color",
                    "name" = "class") |>
      dplyr::arrange(.data[["index"]])
    
  } else {
    file_sf <- file_sf |>
      dplyr::rename("name" = "class") |>
      dplyr::left_join(style |> dplyr::select(dplyr::all_of(c(
        "name", "index"
      )))) |>
      dplyr::mutate(
        class_num = .data[["index"]]
      )
  }
  
  file_bbox <- sf::st_bbox(file_sf)
  
  # Create vector file with `class` converted
  file_gpkg <- fs::file_temp(ext = ".gpkg")
  
  sf::st_write(obj = file_sf, dsn = file_gpkg)
  
  a_srs <- readRDS(url("https://github.com/restore-plus/restore-utils/raw/refs/heads/main/inst/extdata/crs/bdc.rds"))
  
  # Rasterize!
  gdalUtilities::gdal_rasterize(
    src_datasource = file_gpkg,
    dst_filename   = output_file,
    a              = "class_num",
    at             = TRUE,
    tr             = c(res, res),
    te             = file_bbox,
    ot             = "Int16",
    init           = 255,
    a_nodata       = 255,
    co             = c(
      "COMPRESS=ZSTD",
      "PREDICTOR=2",
      "ZSTD_LEVEL=1",
      "BIGTIFF=YES"
    ),
    a_srs          = a_srs
  )
  
  # Save style
  sits:::.colors_qml(color_table = style, file = output_style)
  
  # Return!
  output_file
}

# Step 2.2 -- Get labels' style from ML model
style <- tibble::tibble(
  name = sits_labels(model),
  index = 1:length(sits_labels(model)),
  color = pals::cols25(length(sits_labels(model)))
)

# Step 2.3 -- Aplly rasterize function to all files in directory that has the same version and gpkg extension
to_raster <- paste0("*class_", version, "*.gpkg")
raster_files <- fs::dir_ls(class_dir, glob = to_raster) |>
  purrr::map(function(file) {
    file_name <- fs::path_file(file)
    
    cli::cli_inform("Processing: {file_name}")
    
    sits_rasterize_segments(
      file       = file,
      res        = 10,
      style      = style,
      class_raster_dir = class_raster_output
    )
  })


# ============================================================
# 3. SITS Cube
# ============================================================

# Step 3.1 -- Get labels associated to the trained model data set (Enumerate them in the order they appear according to "sits_labels(model)")
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
  data_dir = class_raster_output, # classified raster file cannot be in the same folder as the classified gpkg file
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
  expected_ua = c("Hidrografia_Rio" = 0.70,
                  "Hidrografia_Lago" = 0.70,
                  "Vegetacao_Natural_Nao_Florestal_Herbacea_Umida" = 0.10,
                  "Vegetacao_Natural_Nao_Florestal_Herbacea_Seca_Mais_Biomassa" = 0.10,
                  "Vegetacao_Natural_Nao_Florestal_Herbacea_Seca_Menos_Biomassa" = 0.10,
                  "Vegetacao_Natural_Nao_Florestal_Vereda" = 0.10,
                  "Vegetacao_Natural_Nao_Florestal_Mata" = 0.10,
                  "Fogo_Recente_Em_Vegetacao_Natural_Nao_Florestal" = 0.10,
                  #"Fogo_Antigo_Em_Vegetacao_Natural_Nao_Florestal" = 0.10,
                  "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Agricultura" = 0.10,
                  #"Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Reservatorio" = 0.10,
                  "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Solo_Exposto" = 0.10,
                  "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Agricultura_Antigo" = 0.10,
                  #"Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Reservatorio_Antigo" = 0.10,
                  "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Solo_Exposto_Antigo" = 0.10
                    ), #comentar a classe caso ela não esteja presente no modelo, senão haverá erro nesta função
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
  alloc = "alloc_50",
  overhead = 1.2, # overproportion to avoid border pixels
  progress = TRUE,
  multicores = 24)

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
samples_validation <- st_read(paste0(samples_dir, "")) #full map validation samples

# Step 5.2 -- Calculate accuracy
area_acc_full_map <- sits_accuracy(cube, 
                          validation = samples_validation,
                          multicores = 28) # adapt to your computer CPU core availability

# Step 5.3 -- Print the area estimated accuracy
area_acc_full_map

# Step 5.4 -- Show confusion matrix
area_acc_full_map$error_matrix

# ============================================================
# 6. Plotting Full Map Accuracy
# ============================================================
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
  scale_y_discrete(limits = rev) +
  #scale_x_discrete(labels = new_label) +
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
                                   mapped_area_ha = round(as.numeric(area_acc_full_map$areas_pixels), 2),
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
mask_label <- c("1" = "Deforestation Mask",
                "0" = "Others")

prodes_mask <- sits_cube(source = "BDC",
                         collection = "SENTINEL-2-16D",
                         data_dir = aux_dir,
                         parse_info = c("tile", "start_date", "end_date", "band", "version"),
                         bands = "class",
                         version = "v2024",
                         labels = mask_label)

cube_reclass <- sits_reclassify(cube = cube,
                                mask = prodes_mask,
                                rules = list("Deforestation" = cube %in% c("Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Agricultura",
                                                                           #"Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Reservatorio",
                                                                           "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Solo_Exposto"
                                                                          ),
                                
                                              "Other Classes" = cube %in% c("Hidrografia_Rio",
                                                                            "Hidrografia_Lago",
                                                                            "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Agricultura_Antigo",
                                                                            #"Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Reservatorio_Antigo",
                                                                            "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Solo_Exposto_Antigo",
                                                                            "Vegetacao_Natural_Nao_Florestal_Herbacea_Umida",
                                                                            "Vegetacao_Natural_Nao_Florestal_Herbacea_Seca_Mais_Biomassa",
                                                                            "Vegetacao_Natural_Nao_Florestal_Herbacea_Seca_Menos_Biomassa",
                                                                            "Vegetacao_Natural_Nao_Florestal_Vereda",
                                                                            "Vegetacao_Natural_Nao_Florestal_Mata",
                                                                            #"Fogo_Antigo_Em_Vegetacao_Natural_Nao_Florestal",
                                                                            "Fogo_Recente_Em_Vegetacao_Natural_Nao_Florestal"
                                                                          )
                                ),
                                multicores = 24,
                                memsize = 180,
                                version = "degradation",
                                output_dir = class_raster_output,
                                progress = TRUE)
plot(cube_reclass)

# 7.2 -- Sampling design prodes
sampling_design <- sits_sampling_design(
  cube = cube_reclass,
  expected_ua = c(
    "Deforestation" = 0.70,
    "Other Classes" = 0.90
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
samples_sf_file_path <- file.path(samples_dir, paste0("samples-validation-desmat_", version, "_", process_version, ".gpkg"))

# 7.6 -- Save samples_sf object as GPKG file
sf::st_write(samples_sf, samples_sf_file_path, append = FALSE)

# ============================================================
# 8. Accuracy assessment of PRODES Adjusted Map classified images
# ============================================================

# Step 8.1 -- Get validation samples points (in geographical coordinates - lat/long)
samples_validation <- st_read(paste0(samples_dir, "")) #prodes adjusted validation samples

# Step 8.2 -- Calculate accuracy
area_acc_prodes <- sits_accuracy(cube_reclass, 
                                 validation = samples_validation,
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
  scale_y_discrete(limits = rev) +
  #scale_x_discrete(labels = new_label) +
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
  scale_y_discrete(limits = rev) +
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
                          mapped_area_ha = round(as.numeric(area_acc_prodes$areas_pixels), 2),
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
  scale_y_discrete(limits = rev) +
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
# 10.1 -- Reclassify classified cube
mask_label <- c("1" = "Deforestation Mask",
                "0" = "Others")

prodes_mask <- sits_cube(source = "BDC",
                         collection = "SENTINEL-2-16D",
                         data_dir = aux_dir,
                         parse_info = c("X1", "X2", "tile", "start_date", "end_date", "band", "version"),
                         bands = "class",
                         version = "v2024",
                         labels = mask_label)




cube_reclass <- sits_reclassify(cube = cube,
                                mask = prodes_mask,
                                rules = list("Deforestation" = cube %in% c("Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Agricultura",
                                                                           #"Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Reservatorio",
                                                                           "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Solo_Exposto"
                                                                          ),
                                             "Degradation" = cube %in% c(#"Fogo_Antigo_Em_Vegetacao_Natural_Nao_Florestal",
                                                                         "Fogo_Recente_Em_Vegetacao_Natural_Nao_Florestal"
                                                                        ),
                                             "Other Classes" = cube %in% c("Hidrografia_Rio",
                                                                           "Hidrografia_Lago",
                                                                           "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Agricultura_Antigo",
                                                                           #"Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Reservatorio_Antigo",
                                                                           "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Solo_Exposto_Antigo",
                                                                           "Vegetacao_Natural_Nao_Florestal_Herbacea_Umida",
                                                                           "Vegetacao_Natural_Nao_Florestal_Herbacea_Seca_Mais_Biomassa",
                                                                           "Vegetacao_Natural_Nao_Florestal_Herbacea_Seca_Menos_Biomassa",
                                                                           "Vegetacao_Natural_Nao_Florestal_Vereda",
                                                                           "Vegetacao_Natural_Nao_Florestal_Mata"
                                                                          )
                                             ),
                                multicores = 24,
                                memsize = 180,
                                version = "degradation",
                                output_dir = class_raster_dir,
                                progress = TRUE)
plot(cube_reclass)

# 10.2 -- Sampling design degradation
sampling_design <- sits_sampling_design(
  cube = cube_reclass,
  expected_ua = c(
    "Deforestation" = 0.70,
    "Degradation" = 0.50, 
    "Other Classes" = 0.95
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
  alloc = "alloc_120",
  overhead = 1.2, # overproportion to avoid border pixels
  progress = TRUE,
  multicores = 28)

# 10.4 -- Total of each class
samples_sf%>% group_by(label) %>% summarise(num = n())

# 10.5 -- Define File Path
samples_sf_file_path <- file.path(samples_dir, paste0("samples-validation-desmat-degrad_", version, "_", process_version, ".gpkg"))

# 10.6 -- Save samples_sf object as GPKG file
sf::st_write(samples_sf, samples_sf_file_path, append = FALSE)

# ============================================================
# 11. Accuracy assessment of PRODES Degradation Adjusted Map classified images
# ============================================================

# Step 11.1 -- Get validation samples points (in geographical coordinates - lat/long)
samples_validation <- st_read(paste0(samples_dir, "")) #prodes adjusted validation samples with degradation classes

# Step 11.2 -- Calculate accuracy
area_acc_prodes <- sits_accuracy(cube_reclass, 
                          validation = samples_validation,
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
  scale_y_discrete(limits = rev) +
  #scale_x_discrete(labels = new_label) +
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
  scale_y_discrete(limits = rev) +
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
                          mapped_area_ha = round(as.numeric(area_acc_prodes$areas_pixels), 2),
                          error_adj_area_ha = round(as.numeric(area_acc_prodes$error_ajusted_area), 2),
                          conf_interval_ha = round(as.numeric(area_acc_prodes$conf_interval), 2)
                         ) %>%
                tidyr::pivot_longer(cols = -class,
                                    names_to = "tipo_area",
                                    values_to = "area")

# Step 10.6 -- Plot Error-Adjusted Area (ha)
ggplot(class_areas, aes(x = tipo_area, y = class, fill = area)) + 
  geom_tile(color = NA, fill = "white") +
  geom_text(aes(label = area), size = 4) +
  geom_vline(xintercept = c(1.5, 2.5), color = "grey30") +
  labs(y = "Class",
       x = "Metrics", 
       title = "Area Metrics") +
  scale_y_discrete(limits = rev) +
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