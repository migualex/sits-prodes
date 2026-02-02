# ============================================================
#  Validation and accuracy measurement
# ============================================================

## I. Load Required Libraries
library(tibble)
library(sits)
library(terra)
library(sf)
library(dplyr)

## II. Define the date and time for the start of processing
date_process <- format(Sys.Date(), "%Y-%m-%d_")
time_process <- format(Sys.time(), "%Hh%Mm", tz = "America/Sao_Paulo")
process_version <- paste0(date_process, time_process)

## III. Define the paths for files and folders needed in the processing
data_dir <- "data/class"
output_dir <- "data/class-raster" # classified raster file cannot be in the same folder as the classified gpkg file
samples_dir <- "data/raw/samples"
model <- readRDS("~/sits-prodes/prodes.amz/data/rds/model/random_forest/RF-model_4-tiles-012015-012014-013015-013014_3y-period-2022-07-28_2025-07-28_with-df-mask-with-all-samples_2026-01-22_09h35m.rds")
version <- "rf-3y-012014-with-df-mask"


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
raster_files <- fs::dir_ls(data_dir, glob = "*class_rf-3y-012014-with-df-mask*.gpkg") |>
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
# 3. Accuracy assessment of classified images
# ============================================================

# Step 3.1 -- Get validation samples points
samples_validation <- st_read("~/sits-prodes/prodes.amz/data/raw/samples/samples_validation_points.gpkg")

# Step 3.2 -- Filter validation samples according to the classified tile
samples_validation_012014 <- samples_validation %>%
  filter(tile == "012014")

# Step 3.3 -- Calculate accuracy
area_acc <- sits_accuracy(cube, 
                          validation = samples_validation_012014,
                          multicores = 2) # adapt to your computer CPU core availability

# Step 3.4 -- Print the area estimated accuracy
area_acc

# Step 3.5 -- Show confusion matrix
area_acc$error_matrix