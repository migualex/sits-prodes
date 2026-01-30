# ============================================================
#  Validation and accuracy measurement
# ============================================================

## I. Load Required Libraries
library(tibble)
library(sits)
library(terra)
library(sf)
library(dplyr)
library(mapview)

## II. Define the date and time for the start of processing
date_process <- format(Sys.Date(), "%Y-%m-%d_")
time_process <- format(Sys.time(), "%Hh%Mm", tz = "America/Sao_Paulo")
process_version <- paste0(date_process, time_process)

## III. Define the paths for files and folders needed in the processing
data_dir <- "data/class"
output_dir <- "data/class-raster" # classified raster file cannot be in the same folder as the classified gpkg file
samples_dir <- "data/raw/samples"
model <- readRDS("~/sits-prodes/prodes.amz/data/rds/model/random_forest/2025_12_12_08h46m_RF_1y_012015_012014_013015_013014_with-df-mask.rds")
version <- "rf-1y-013015-with-df-mask"


# ============================================================
# 1. Create raster file from classified vector map
# ============================================================

# 1.1 -- Function to rasterize
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

# 1.2 -- Get labels' style from ML model
style <- tibble::tibble(
  name = sits_labels(model),
  index = 1:length(sits_labels(model)),
  color = pals::cols25(length(sits_labels(model)))
)

# 1.3 -- Aplly rasterize function to all files in directory that has the same version and gpkg extension
raster_files <- fs::dir_ls(data_dir, glob = "*class_rf-1y-013015-with-df-mask*.gpkg") |>
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
# 2. SITS Cube
# ============================================================

# Step 2.1 -- Get labels associated to the trained model data set
sits_labels(model)

# Step 2.2 -- Labels of the classified image (Enumerate them in the order they appear according to "sits_labels(model)")
labels <- c(
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

# Step 2.3 -- Load the original cube with classified raster file
cube <- sits_cube(
  source = "BDC",
  collection = "SENTINEL-2-16D",
  bands = "class",
  labels = labels,
  data_dir = output_dir, # classified raster file cannot be in the same folder as the classified gpkg file
  version = version,
  parse_info = c("satellite", "sensor", "tile", "start_date", "end_date", 
                 "band", "version"))

# Step 2.4 -- Extract and define some information
tiles_class <- paste(cube$tile, collapse = "-")
dates <- sits_timeline(cube)
no.years <- paste0(floor(lubridate::interval(dates[1], dates[length(dates)]) / lubridate::years(1)), "y")
var <- "with-df-mask"
version <- paste("rf", no.years, tiles_class, var, sep = "-")

# ============================================================
# 3. Stratified random sampling
# ============================================================

# 3.1 -- Sampling design
sampling_design <- sits_sampling_design(
  cube = cube,
  expected_ua = c(
    "AGUA" = 0.95,
    "DESMAT_ARVORE_REMANESCE" = 0.40, 
    "DESMAT_CORTE_RASO" = 0.75, 
    "DESMAT_CORTE_RASO_DM" = 0.9,  
    "DESMAT_DEGRAD_FOGO" = 0.70, 
    "DESMAT_VEG" = 0.70,  
    "DESMAT_VEG_DM" = 0.75, 
    "FLO_DEGRAD" = 0.70, 
    "FLO_DEGRAD_FOGO" = 0.70,
    "FLORESTA" = 0.95,
    "NF" = 0.90,
    "ROCHA" = 0.80,
    "WETLANDS" = 0.70
  ),
  alloc_options = c(120, 100, 75, 50, 30),
  std_err = 0.01,
  rare_class_prop = 0.01
)

# 3.2 -- Show sampling design
sampling_design

# 3.3 -- Generate stratified samples
samples_sf <- sits_stratified_sampling(
  cube = cube,
  sampling_design = sampling_design,
  alloc = "alloc_50",
  overhead = 1.2, #proporção de pixels a mais que eu quero que pegue para descatar pixels de borda
  progress = TRUE,
  multicores = 12)

samples_sf%>% group_by(label) %>% summarise(num = n())


# 3.4 -- Define File Path
samples_sf_file_path <- file.path(samples_dir, paste0("samples-validation_", version, "_", process_version, ".gpkg"))

# 3.4 -- Save samples_sf object as GPKG file
sf::st_write(samples_sf, samples_sf_file_path, append = FALSE)

# ============================================================
# 4. Accuracy assessment of classified images
# ============================================================

# Step 4.1 -- Get validation samples points (coordenadas geograficas)
samples_validation <- st_read(samples_sf_file_path)

# Step 4.2 -- Calculate accuracy
area_acc <- sits_accuracy(cube, 
                          validation = samples_validation,
                          multicores = 10) # adapt to your computer CPU core availability

# Step 4.3 -- Print the area estimated accuracy
area_acc

# Step 4.4 -- Show confusion matrix
area_acc$error_matrix