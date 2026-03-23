# ============================================================
#  Stratification of validation samples
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
aux_dir          <- "data/raw/auxiliary"
version          <- "rf-1y-013014-all-classes"

# Step 1.4 -- Create the directory for storing class rasters, including any necessary parent directories. Suppress warnings if the directory already exists.
dir.create(class_raster_dir, recursive = TRUE, showWarnings = FALSE)

# Step 1.5 -- Get the list of validation sample files matching the version pattern in the samples directory
samples_validation_list <- dir(
  samples_dir,
  pattern = paste0(".*", version, ".*\\.gpkg$"),
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
                 "band", "version"))

# ============================================================
# 4. Full Map Stratified Random Sampling
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

# 5.3 -- Extract the cube REQUIRED
cube_reclass <- cube_reclass[[1]]

# 5.4 -- Sampling design degradation
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

# 5.5 -- Show sampling design
sampling_design

# 5.6 -- Generate stratified samples
samples_sf <- sits_stratified_sampling(
  cube = cube_reclass,
  sampling_design = sampling_design,
  alloc = "alloc_100",
  overhead = 1.2, # overproportion to avoid border pixels
  progress = TRUE,
  multicores = 24)

# 5.7 -- Total of each class
samples_sf%>% group_by(label) %>% summarise(num = n())

# 5.8 -- Define File Path
samples_sf_file_path <- file.path(samples_dir, paste0("samples-validation-desmat-degrad_", version, "_", process_version, ".gpkg"))

# 5.9 -- Save samples_sf object as GPKG file
sf::st_write(samples_sf, samples_sf_file_path, delete_dsn = TRUE, append = FALSE)