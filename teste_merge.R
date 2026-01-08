# Packages
library(sits)
library(torch)
library(sf)
library(dplyr)
library(glue)

# Parameters
output_dir <- "data/"
roi_name   <- "SC"

# Cube dates
start_date <- "2024-07-01"
end_date   <- "2025-08-12"

# Cube bands
bands <- c("VV","VH")

# Regularization period
period <- "P16D"

# Hardware
memsize <- 10
multicores <- 1

# Get ROI
#roi_sf <- geobr::read_state(roi_name)

roi <- c(
  lon_min = -64.20726838394573,
  lat_min = -8.55912145745303,
  lon_max = -63.20391740487580,
  lat_max = -7.58418972216315
)

cube <- sits_cube(
  source     = "MPC",
  collection = "SENTINEL-1-RTC",
  bands      = bands,
  start_date = start_date,
  end_date   = end_date,
  roi        = roi
)

# Get valid timeline (Internal SITS function)
timeline <- sits:::.gc_get_valid_timeline(cube, period)

# Define output dir
output_dir_s2 <- fs::path(output_dir) / "s2-reg"

# Create output dir
fs::dir_create(output_dir_s2)

# Regularize tile by tile
slider::slide_dfr(cube, function(row) {
  # Wait few seconds before move to the next tile
  Sys.sleep(10)
  
  # Inform user about the current tile
  print(glue::glue("Processing: {row[['tile']]}"))
  
  # Process!
  tryCatch(
    {
      sits_regularize(
        cube       = row,
        period     = period,
        res        = 10,
        roi        = roi_sf,
        multicores = memsize,
        output_dir = output_dir_s2,
        timeline   = timeline
      )
    },
    error = function(e) { NULL }
  )
})

# Reload cube
cube_local <- sits_cube(
  source     = "MPC",
  collection = "SENTINEL-2-L2A",
  data_dir   = output_dir_s2
)

# Check cube bands (must be "TRUE")
all(sits_bands(cube_local) %in% setdiff(bands, "CLOUD"))

# Check timeline
sits_timeline(cube_local)

# Check tiles
nrow(cube_local)

# Check if cube is regular (must be "TRUE")
sits:::.cube_is_regular(cube_local)