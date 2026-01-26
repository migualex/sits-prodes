# ============================================================
#  Validation and accuracy measurement
# ============================================================

## I. Load Required Libraries
library(tibble)
library(sits)
library(terra)
library(sf)

## II. Define the date and time for the start of processing
date_process <- format(Sys.Date(), "/%Y_%m_%d_")
time_process <- format(Sys.time(), "%Hh%Mm_", tz = "America/Sao_Paulo")
process_version <- paste0(date_process, time_process)

## III. Define the paths for files and folders needed in the processing
data_dir <- "data/class"

#val_samples_dir <- "data/raw/validation_samples"
#train_samples_dir <- "data/raw/training_samples"


# ============================================================
# 1. ML Model and SITS Cube
# ============================================================

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
  "12" = "WETLANDS"
)

# Step 1.1 -- Load the original cube
cube <- sits_cube(
  source = "BDC",
  collection = "SENTINEL-2-16D",
  bands = "class",
  labels = labels,
  data_dir = data_dir,
  version = "rf-3y-012014-with-df-mask",
  parse_info = c("satellite", "sensor", "tile", "start_date", "end_date", 
                 "band", "version"))

# Step 1.2 -- Extract and define some information
tiles_class <- paste(cube$tile, collapse = "-")
dates <- sits_timeline(cube)
no.years <- paste0(floor(lubridate::interval(dates[1], dates[length(dates)]) / lubridate::years(1)), "y")
var <- "with-df-mask"
version <- paste("rf", no.years, tiles_class, var, sep = "-")

# Step 1.3 -- Recovery ML Model
rf_model <- readRDS("data/rds/model/random_forest/RF-model_4-tiles-012015-012014-013015-013014_3y-period-2022-07-28_2025-07-28_with-df-mask-with-all-samples_2026-01-22_09h35m.rds")


# ============================================================
# 2. Classification Map
# ============================================================

# Step 2.1 -- Read GeoPackage File
class_data <- st_read(classified_gpkg_path)

# Step 2.2 -- Verify unique classes
unique(class_data$class)

# Step 2.3 -- Define the classes of the probability cube
class_mapping <- c(
  "DESMAT_CORTE_RASO" = 1,
  "DESMAT_CORTE_RASO_DM" = 2,
  "DESMAT_ARVORE_REMANESCE" = 3,
  "FLORESTA" = 4,
  "FLO_DEGRAD_FOGO" = 5,
  "FLO_DEGRAD" = 6,
  "DESMAT_VEG" = 7,
  "DESMAT_VEG_DM" = 8,
  "AGUA" = 9,
  "NF" = 10,
  "ROCHA" = 11,
  "WETLANDS" = 12
)

# Step 2.4 -- Add numeric column
class_data$class_num <- class_mapping[class_data$class]

# Step 2.5 -- Verify if all of them were mapped
table(class_data$class, class_data$class_num, useNA = "ifany")

# ============================================================
# 2. Create raster
# ============================================================

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

#map_class <- read_sf("data/class/SENTINEL-2_MSI_012014_2022-07-28_2025-07-28_class_rf-3y-012014-with-df-mask.gpkg")

model <- readRDS("~/sits-prodes/prodes.amz/data/rds/model/random_forest/RF-model_4-tiles-012015-012014-013015-013014_3y-period-2022-07-28_2025-07-28_with-df-mask-with-all-samples_2026-01-22_09h35m.rds")

style <- tibble::tibble(
  name = sits_labels(model),
  index = 1:length(sits_labels(model)),
  color = pals::cols25(length(sits_labels(model)))
)

input_dir <- "data/class/"
output_dir <- "data/class-raster"

raster_files <- fs::dir_ls(input_dir, glob = "*class_rf-3y-012014-with-df-mask*.gpkg") |>
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
# 3. 
# ============================================================

# Step 3.1 -- Sampling design
sampling_design <- sits_sampling_design(
  cube = class_datacube,
  expected_ua = c(
    "DESMAT_CORTE_RASO" = 0.75,
    "DESMAT_CORTE_RASO_DM" = 0.70,
    "DESMAT_ARVORE_REMANESCE" = 0.70, 
    "DESMAT_VEG_DM" = 0.70, 
    "FLO_DEGRAD_FOGO" = 0.70,
    "FLO_DEGRAD" = 0.70,
    "NF" = 0.70,
    "ROCHA" = 0.70,
    "FLORESTA" = 0.75,  
    "DESMAT_VEG" = 0.70,  
    "AGUA" = 0.70, 
    "WETLANDS" = 0.70
  ),
  alloc_options = c(120, 100),
  std_err = 0.01,
  rare_class_prop = 0.1
)


# Step 3.2 -- Show sampling design
sampling_design

# Step 3.3 -- Generate stratified random samples
samples_sf <- sits_stratified_sampling(
  cube = class_datacube,
  sampling_design = sampling_design,
  alloc = "alloc_120",
  multicores = 4 # adapt to your computer CPU core availability
  )

# Step 3.4 -- Save samples in a shapefile
sf::st_write(samples_sf, 
             file.path(val_samples_dir, "validation_samples.shp"), 
             append = FALSE
             )


# ============================================================
# 4. Accuracy assessment of classified images
# ============================================================

# Step 4.5 -- Get ground truth points
ground_truth <- system.file(
  "class/samples_validation.csv", package = "sitsdata"
  )

# Step 2.6 -- Calculate accuracy according to Olofsson's method
area_acc <- sits_accuracy(probs_datacube_class, 
                          validation = ground_truth,
                          multicores = 4 # adapt to your computer CPU core availability
                          )





model <- readRDS("~/sits-prodes/prodes.amz/data/rds/model/random_forest/RF-model_4-tiles-012015-012014-013015-013014_3y-period-2022-07-28_2025-07-28_with-df-mask-with-all-samples_2026-01-22_09h35m.rds")

version <- "rf-3y-012014-with-df-mask"

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

class_cube <- sits_cube(
  source = "BDC",
  collection = "SENTINEL-2-16D",
  bands = "class",
  labels = labels,
  data_dir = data_dir,
  version = version,
  parse_info = c("satellite", "sensor", "tile", "start_date", "end_date", 
                 "band", "version")
)


validation <- readRDS("data/rds/time_series/samples-val_4-tiles-012015-012014-013015-013014_3y-period-2022-07-28_2025-07-28_with-df-mask_2026-01-21_10h10m.rds")


validation_csv <- validation |>
  dplyr::select(latitude, longitude, label)

write.csv(
  validation_csv,
  "data/test/validation_accuracy.csv",
  row.names = FALSE
)


# Calculate accuracy according to Olofsson's method
area_acc <- sits_accuracy(class_cube, 
                          validation = validation_csv,
                          multicores = 4)

# Step 2.7 -- Print the area estimated accuracy
area_acc


# Step 2.8 -- Show confusion matrix
area_acc$error_matrix


validation_clean <- validation |>
  dplyr::select(
    -time_series,
    -cube
  )

