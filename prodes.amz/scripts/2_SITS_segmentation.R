# ============================================================
#  Segmentation from Linear Spectral Mixture Model (LSMM) for features training and classification
# ============================================================

## I. Load Required Libraries
library(sits)

## II. Define the paths for files and folders needed in the processing
mixture_path <- "data/raw/mixture_model"
segments_path <- "data/segments"

date_process <- format(Sys.Date(), "%Y-%m-%d")

# ============================================================
# 1. Load Mixture Model Data Cubes from a collection
# ============================================================

mm_cube_local <- sits_cube(
  source = "BDC",
  tiles = '012014',
  collection = "SENTINEL-2-16D",
  data_dir = mixture_path,
  bands = c("SOIL", "VEG", "WATER"),
  start_date = '2024-08-01',
  end_date = '2025-07-31',
  progress = TRUE)

# ============================================================
# 2. Segment the mixture model cube based on Fraction Images features
# ============================================================

version <- paste("comp05-hex", date_process)

sits_segment_start <- Sys.time()
mm_cube_segments <- sits_segment(
  cube = mm_cube_fraction_features,
  seg_fn = sits_snic(
    grid_seeding = "hexagonal",
    spacing = 10,
    compactness = 0.5,
    padding = 0
  ),
  memsize = 80, #adapt to your computer memory availability
  multicores = 4, #adapt to your computer CPU core availability
  output_dir = segments_path,
  version = version
)
sits_segment_end <- Sys.time()
process_duration_sits_segment <- round(sits_segment_end-sits_segment_start,2)
process_duration_sits_segment

print("Segmentation has finished")