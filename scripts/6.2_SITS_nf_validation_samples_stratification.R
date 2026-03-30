# ============================================================
#  Stratification of Non-Forest validation samples
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
aux_dir          <- "data/raw/auxiliary"
version          <- "rf-2y-novos-segmentos"

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
to_raster <- paste0(".*_class_rf-2y-", "01[4-6]00[1-3]", "-novos-segmentos*\\.gpkg$")

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
  tile_period_dir <- file.path(
    class_dir,
    tile_id,
    "raster"
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
cube_dirs <- list.dirs(class_dir, recursive = TRUE)

cube_dirs <- cube_dirs[
  sapply(cube_dirs, function(x) {
    files <- list.files(x, pattern = "\\.tif$")
    any(grepl(version, files))
  })
]

cube_dirs <- grep("*raster", cube_dirs, value = TRUE)

labels <- c(
  x = sits_labels(model)
)
names(labels) <- 1:length(labels)

# Step 3.2 -- Load the original cube with classified raster file
cube <- sits_cube(
  source     = "BDC",
  collection = "SENTINEL-2-16D",
  bands      = "class",
  labels     = labels,
  data_dir   = cube_dirs, # classified raster file cannot be in the same folder as the classified gpkg file
  version    = version,
  parse_info = c("satellite", "sensor", "tile", "start_date", "end_date",
                 "band", "version"))

# ============================================================
# 4. Full Map Stratified Random Sampling
# ============================================================

# 4.1 -- Sampling design
sampling_design <- sits_sampling_design(
  cube = cube,
  expected_ua = c(
    "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Agricultura" = 0.70,
    "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Solo_Exposto" = 0.70,
    "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Agricultura_Antigo" = 0.70,
    "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Solo_Exposto_Antigo" = 0.70,
    "Hidrografia_Lago" = 0.95,
    "Hidrografia_Rio" = 0.95,
    "Fogo_Recente_Em_Vegetacao_Natural_Nao_Florestal" = 0.90,
    "Vegetacao_Natural_Nao_Florestal_Herbacea_Seca_Mais_Biomassa" = 0.90,
    "Vegetacao_Natural_Nao_Florestal_Herbacea_Seca_Menos_Biomassa" = 0.90,
    "Vegetacao_Natural_Nao_Florestal_Herbacea_Seca_Pos_Fogo" = 0.90,
    "Vegetacao_Natural_Nao_Florestal_Herbacea_Umida" = 0.90,
    "Vegetacao_Natural_Nao_Florestal_Transicao_Florestal" = 0.90,
    "Vegetacao_Natural_Nao_Florestal_Mata" = 0.90,
    "Vegetacao_Natural_Nao_Florestal_Vereda" = 0.90
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
                         version = "lavrado-v2024",
                         labels = mask_label)

cube_reclass <- sits_reclassify(
  cube = cube,
  mask = prodes_mask,
  rules = list(
    "Deforestation" =
      cube %in% c(
        "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Agricultura",
        "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Solo_Exposto",
        "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Agricultura_Antigo",
        "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Solo_Exposto_Antigo"
      ),
    "Water" =
      cube %in% c(
        "Hidrografia_Lago",
        "Hidrografia_Rio"
      ),
    "Grassland" =
      cube %in% c(
        "Fogo_Recente_Em_Vegetacao_Natural_Nao_Florestal",
        "Vegetacao_Natural_Nao_Florestal_Herbacea_Seca_Mais_Biomassa",
        "Vegetacao_Natural_Nao_Florestal_Herbacea_Seca_Menos_Biomassa",
        "Vegetacao_Natural_Nao_Florestal_Herbacea_Seca_Pos_Fogo",
        "Vegetacao_Natural_Nao_Florestal_Herbacea_Umida",
        "Vegetacao_Natural_Nao_Florestal_Transicao_Florestal"
      ),
    "Forest" = 
      cube %in% c(
        "Vegetacao_Natural_Nao_Florestal_Mata",
        "Vegetacao_Natural_Nao_Florestal_Vereda"
      )
  ),
  multicores = 24,
  memsize = 180,
  version = paste("prodes-degradation", version, sep = "-"),
  output_dir = "data/class/all_tiles",
  progress = TRUE
)

# 5.4 -- Sampling design degradation
sampling_design <- sits_sampling_design(
  cube = cube_reclass,
  expected_ua = c(
    "Deforestation" = 0.70,
    "Water" = 0.95, 
    "Grassland" = 0.90,
    "Forest" = 0.90
  ),
  alloc_options = c(120, 100, 75, 50, 30),
  std_err = 0.01,
  rare_class_prop = 0.04
)

# 5.5 -- Show sampling design
sampling_design

# 5.6 -- Generate stratified samples
samples_sf <- sits_stratified_sampling(
  cube = cube_reclass,
  sampling_design = sampling_design,
  alloc = "alloc_120",
  overhead = 1.2, # overproportion to avoid border pixels
  progress = TRUE,
  multicores = 24)

# 5.7 -- Total of each class
samples_sf%>% group_by(label) %>% summarise(num = n())

# 5.8 -- Define File Path
samples_sf_file_path <- file.path(samples_dir, paste0("samples-validation-desmat-degrad_", version, "_", process_version, ".gpkg"))

# 5.9 -- Save samples_sf object as GPKG file
sf::st_write(samples_sf, samples_sf_file_path, delete_dsn = TRUE, append = FALSE)
