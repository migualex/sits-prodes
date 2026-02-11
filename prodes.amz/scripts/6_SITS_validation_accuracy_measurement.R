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

pastas <- c("data/class/raster", "data/raw/aux")
lapply(pastas, dir.create, recursive = FALSE, showWarnings = FALSE)

# Step 1.3 -- Define the paths for files and folders needed in the processing
data_dir <- "data/class"
output_dir <- "data/class/raster" # classified raster file cannot be in the same folder as the classified gpkg file
samples_dir <- "data/raw/samples"
plots_path    <- "data/plots/"
aux_dir <- "data/raw/aux"
#model <- readRDS("~/sits-prodes/prodes.amz/data/rds/model/random_forest/RF-model_4-tiles-012015-012014-013015-013014_1y-period-2024-07-27_2025-07-28_no-remaining-trees_2026-02-02_15h14m.rds")
version <- "rf-1y-012014-all-classes"


# ============================================================
# 2. Create raster file from classified vector map
# ============================================================

# Step 2.1 -- Function to rasterize
sits_rasterize_segments <- function(file, res, output_dir, style = NULL) {
  stopifnot(!is.null(res))
  stopifnot(!is.null(file))
  stopifnot(!is.null(output_dir))
  
  # create output dir
  fs::dir_create(output_dir, recurse = TRUE)
  
  # expand paths
  file <- fs::path_expand(file)
  output_dir <- fs::path_expand(output_dir)
  
  # define output files
  output_file_base <- fs::path(output_dir) / fs::path_file(file) |>
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
raster_files <- fs::dir_ls(data_dir, glob = "*class_rf-1y-012014-no-remaining-trees*.gpkg") |>
  purrr::map(function(file) {
    file_name <- fs::path_file(file)
    
    cli::cli_inform("Processing: {file_name}")
    
    sits_rasterize_segments(
      file       = file,
      res        = 10,
      style      = style,
      output_dir = output_dir
    )
  })


# ============================================================
# 3. SITS Cube
# ============================================================

# Step 3.1 -- Get labels associated to the trained model data set
sits_labels(model)

# Step 3.2 -- Labels of the classified image (Enumerate them in the order they appear according to "sits_labels(model)")
labels <- c(
#   "1" = "AGUA",
#   "2" = "DESMAT_CORTE_RASO",
#   "3" = "DESMAT_CORTE_RASO_DM",
#   "4" = "DESMAT_DEGRAD_FOGO",
#   "5" = "DESMAT_VEG",
#   "6" = "DESMAT_VEG_DM",
#   "7" = "FLO_DEGRAD", 
#   "8" = "FLO_DEGRAD_FOGO",
#   "9" = "FLORESTA",
#   "10" = "NF",
#   "11" = "ROCHA",
#   "12" = "WETLANDS"
# )  
"1" = "AGUA",
"2" = "DESMAT_ARVORE_REMANESCE",
"3" = "DESMAT_CORTE_RASO",
"4" = "DESMAT_CORTE_RASO_DM",
"5" = "DESMAT_DEGRAD_FOGO",
"6" = "DESMAT_VEG",
"7" = "DESMAT_VEG_DM",
"8" = "FLO_DEGRAD",
"9" = "FLO_DEGRAD_FOGO",
"10" = "FLORESTA",
"11" = "NF",
"12" = "ROCHA",
"13" = "WETLANDS"
)

# Step 3.3 -- Load the original cube with classified raster file
cube <- sits_cube(
  source = "BDC",
  collection = "SENTINEL-2-16D",
  bands = "class",
  labels = labels,
  data_dir = output_dir, # classified raster file cannot be in the same folder as the classified gpkg file
  version = version,
  parse_info = c("satellite", "sensor", "tile", "start_date", "end_date", 
                 "band", "version")
  )

# Step 3.4 -- Extract and define some information
tiles <- paste(cube$tile, collapse = "-")
dates <- sits_timeline(cube)
no.years <- paste0(floor(lubridate::interval(dates[1], dates[length(dates)]) / lubridate::years(1)), "y")
var <- "no-remaining-trees"
version <- paste("rf", no.years, tiles, var, sep = "-")


# ============================================================
# 4. Full Map Stratified random sampling
# ============================================================

# 4.1 -- Sampling design
sampling_design <- sits_sampling_design(
  cube = cube,
  expected_ua = c(
    "AGUA" = 0.95,
    "DESMAT_ARVORE_REMANESCE" = 0.10, 
    "DESMAT_CORTE_RASO" = 0.70, 
    "DESMAT_CORTE_RASO_DM" = 0.85, # ok
    "DESMAT_DEGRAD_FOGO" = 0.70, 
    "DESMAT_VEG" = 0.70,
    "DESMAT_VEG_DM" = 0.85, # ok
    "FLO_DEGRAD" = 0.70, 
    "FLO_DEGRAD_FOGO" = 0.70,
    "FLORESTA" = 0.95,
    "NF" = 0.70,
    "ROCHA" = 0.10,
    "WETLANDS" = 0.70
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
samples_sf_file_path <- file.path(samples_dir, paste0("samples-validation_", version, "_", process_version, ".gpkg"))

# 4.6 -- Save samples_sf object as GPKG file
sf::st_write(samples_sf, samples_sf_file_path, append = FALSE)

# ============================================================
# 5. Accuracy assessment of Full Map classified images
# ============================================================

# Step 5.1 -- Get validation samples points (in geographical coordinates - lat/long)
samples_validation <- st_read("data/raw/samples_validation/samples-validation_012014_no-remaining_trees_1y_27-07-2024_12-07-2025.gpkg")

# Step 5.2 -- Calculate accuracy
area_acc <- sits_accuracy(cube, 
                          validation = samples_validation,
                          multicores = 28) # adapt to your computer CPU core availability

# Step 5.3 -- Print the area estimated accuracy
area_acc

# Step 5.4 -- Show confusion matrix
area_acc$error_matrix

# ============================================================
# 6. Plotting Full Map Accuracy
# ============================================================

# Step 6.1 -- Change labels' names for plotting purpose
new_label <-c("AGUA" = "Water",
              "DESMAT_ARVORE_REMANESCE" = "Clear Cut Trees", 
              "DESMAT_CORTE_RASO" = "Clear Cut Bar Soil", 
              "DESMAT_CORTE_RASO_DM" = "Old Clear Cut Bare Soil",
              "DESMAT_DEGRAD_FOGO" = "Fire Suppression", 
              "DESMAT_VEG" = "Clear Cut Vegetation",
              "DESMAT_VEG_DM" = "Old Clear Cut Vegetation",
              "FLO_DEGRAD" = "Degradation", 
              "FLO_DEGRAD_FOGO" = "Fire Degradation",
              "FLORESTA" = "Forest",
              "NF" = "Non Forest Natural Vegetation",
              "ROCHA" = "Rock",
              "WETLANDS" = "Wetland"
              )

# Step 6.2 -- Create a tibble from error matrix
matriz_conf <- tibble(as.data.frame(area_acc$error_matrix))

# Step 6.3 -- Plot error matrix
ggplot(matriz_conf, aes(x = Var2, y = Var1, fill = Freq)) +
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
  filename = paste(version, "matrix-confusion.png", sep = "_"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# Step 6.4 -- Convert accuracies results to a data frame
acuracias <- data.frame(class = names(area_acc$accuracy[[1]]),
                        user_accuracy = round(as.numeric(area_acc$accuracy[[1]]), 2),
                        prod_accuracy = round(as.numeric(area_acc$accuracy[[2]]), 2)
                        ) %>% 
                        tidyr::pivot_longer(cols = -class,
                                            names_to = "tipo_acuracia",
                                            values_to = "acuracia")

# Step 6.5 -- Plot accuracies metrics
ggplot(acuracias, aes(x = tipo_acuracia, y = class, fill = acuracia)) + 
  geom_tile(color = "white") +
  geom_text(aes(label = acuracia), size = 4) +
  scale_fill_distiller(palette = "RdYlGn",
                       direction = 1,
                       name = "Value") +
  labs(y = "Class",
       x = "Accuracy", 
       title = "Accuracies",
       caption = paste0("Global Accuracy: ", round(area_acc$accuracy[[3]], 2))) +
  scale_y_discrete(limits = rev,
                   labels = new_label) +
  scale_x_discrete(labels = c("prod_accuracy" = "Prod Acc",
                              "user_accuracy" = "User Acc")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 1, size = 11, margin = margin(t = 10)),
        panel.grid = element_blank())

ggsave(
  filename = paste(version, "accuracies.png", sep = "_"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# Step 6.6 -- Plot Error-Adjusted Area (ha)
class_areas <- data.frame(class = names(area_acc$area_pixels),
                          mapped_area_ha = round(as.numeric(area_acc$areas_pixels), 2),
                          error_adj_area_ha = round(as.numeric(area_acc$error_ajusted_area), 2),
                          conf_interval_ha = round(as.numeric(area_acc$conf_interval), 2)
                          ) %>%
                          tidyr::pivot_longer(cols = -class,
                                              names_to = "tipo_area",
                                              values_to = "area")

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
  filename = paste(version, "areas.png", sep = "_"),
  path = plots_path,
  scale = 1,
  width = 3529,
  height = 1578,
  units = "px",
  dpi = 350,
)

# ============================================================
# 6. PRODES Adjusted Map Accuracy
# ============================================================
# 4.1 -- Reclassify classified cube
mask_label <- c("1" = "Deforestation Mask",
                "0" = "bla")

prodes_mask <- sits_cube(source = "BDC",
                         collection = "SENTINEL-2-16D",
                         data_dir = aux_dir,
                         parse_info = c("X1", "X2", "tile", "start_date", "end_date", "band", "version"),
                         bands = "class",
                         version = "v2024-f",
                         labels = mask_label)

cube_reclass <- sits_reclassify(cube = cube,
                                mask = prodes_mask,
                                rules = list("Deforestation" = cube %in% c("DESMAT_ARVORE_REMANESCE",
                                                                           "DESMAT_CORTE_RASO",
                                                                           "DESMAT_VEG",
                                                                           "DESMAT_DEGRAD_FOGO"
                                                                          ),
                                             "Degradation" = cube %in% c("FLO_DEGRAD",
                                                                         "FLO_DEGRAD_FOGO"
                                                                        ),
                                             "Other Classes" = cube %in% c("AGUA",
                                                                           "DESMAT_CORTE_RASO_DM",
                                                                           "DESMAT_VEG_DM",
                                                                           "FLORESTA",
                                                                           "NF",
                                                                           "ROCHA",
                                                                           "WETLANDS"
                                                                          )
                                             ),
                                multicores = 24,
                                memsize = 180,
                                version = "V4",
                                output_dir = output_dir,
                                progress = TRUE)
plot(cube_reclass)

# 4.2 -- Sampling design degradation
sampling_design <- sits_sampling_design(
  cube = cube_reclass,
  expected_ua = c(
    "Deforestation" = 0.70,
    "Degradation" = 0.50, 
    "Other Classes" = 0.90
  ),
  alloc_options = c(120, 100),
  std_err = 0.01,
  rare_class_prop = 0.05
)

# 4.2 -- Show sampling design
sampling_design

# 4.3 -- Generate stratified samples
samples_sf <- sits_stratified_sampling(
  cube = cube_reclass,
  sampling_design = sampling_design,
  alloc = "alloc_30",
  overhead = 1.2, # overproportion to avoid border pixels
  progress = TRUE,
  multicores = 12)

# 4.4 -- Total of each class
samples_sf%>% group_by(label) %>% summarise(num = n())

# 4.5 -- Define File Path
samples_sf_file_path <- file.path(samples_dir, paste0("samples-validation_", version, "_", process_version, ".gpkg"))

# 4.6 -- Save samples_sf object as GPKG file
sf::st_write(samples_sf, samples_sf_file_path, append = FALSE)

# ============================================================
# 5. Accuracy assessment of Full Map classified images
# ============================================================

# Step 5.1 -- Get validation samples points (in geographical coordinates - lat/long)
samples_validation <- st_read("data/raw/samples_validation/samples-validation_012014_no-remaining_trees_1y_27-07-2024_12-07-2025.gpkg")

# Step 5.2 -- Calculate accuracy
area_acc <- sits_accuracy(cube_reclass, 
                          validation = samples_validation,
                          multicores = 28) # adapt to your computer CPU core availability

# Step 5.3 -- Print the area estimated accuracy
area_acc

# Step 5.4 -- Show confusion matrix
area_acc$error_matrix