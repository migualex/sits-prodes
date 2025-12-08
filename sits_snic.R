# Script to segment Sentinel-2 data cube 012014 from BDC

library(sits)

# path to cube images:
images_path <- "~/PRODES-SITS/cube"
mixture_path <- "~/PRODES-SITS/mixture"
# path to reduced cubes
segments_path <- "~/PRODES-SITS/segments"

# ---- Step 1. Build a local copy of tile cube 012014 from BDC
# access to BDC
cube <- sits_cube(
    source = "BDC",
    collection = "SENTINEL-2-16D",
    bands = c('B02', 'B03', 'B04', 'B05', 'B06','B07', 'B08','B8A', 'B11','B12', 'NDVI', 'NBR','EVI'),
    tiles = c('012014'),  # '012015', '013014', '013015'),
    start_date = '2024-07-01',
    end_date = '2025-08-12',
    progress = TRUE)

# Copy BDC cube to local files
cube_local <- sits_cube_copy(
    cube = cube,
    multicores = 2,
    output_dir = images_path,
    progress = TRUE
)
# Retrieve data cube from local files
cube_012014_local <- sits_cube(
    source = "BDC",
    collection = "SENTINEL-2-16D",
    data_dir = images_path
)

# ----- Part 2. Generation of Random Forest model
# retrieval of time series based on samples by PRODES team
samples_shp <- "~/PRODES-SITS/samples/amostras_segment_mlme291025.shp"

samples <- sits_get_data(
    cube = cube_012014_local,
    samples = samples_shp,
    n_samples_pol = 16,
    multicores = 4,
    progress = TRUE
)
# Build an RF model based on time series samples
rf_model <- sits_train(
    samples = samples,
    ml_method = sits_rfor()
)
# save RF model
saveRDS(rf_model, file = "~/PRODES-SITS/models/rf_model_13bands.rds")

# -- Part 3. Generation of mixture model cube
# mixture model values
end_members <- tibble::tribble(
    ~class,     ~B03,    ~B04,    ~B08,    ~B11,
    "veg",       620,     207,    4450,    3123,
    "soil",     1275,    1501,    3958,    4035,
    "water",     329,      72,      86,      72
)
# create a mixture model cube
mm_cube <- sits_mixture_model(
    data = cube_012014_local,
    endmembers = end_members,
    multicores = 3,
    memsize = 12,
    output_dir = mixture_path
)
# retrieve the local cube with end members from
# mixture model
mm_cube_local <- sits_cube(
    source = "BDC",
    collection = "SENTINEL-2-16D",
    data_dir = mixture_path
)
# select only the end members
mm_cube_end_members <- sits_select(
    data = mm_cube_local,
    bands = c("SOIL", "VEG", "WATER")
)

# --- Part 4. Segmentation of end members cube

# create segments form end members cube
mm_cube_segments <- sits_segment(
    cube = mm_cube_end_members,
    seg_fn = sits_snic(
        grid_seeding = "rectangular",
        spacing = 30,
        compactness = 0.5,
        padding = 15
    ),
    memsize = 12,
    multicores = 4,
    output_dir = segments_path,
    version = "snic-30-05"
)

# --- Part 5. Classification of segments

# Build a vector cube combining segments produced by mixture model
# with local data cube of bands and indices (10 bands + 3 indices)
local_segs_cube <- sits_cube(
    source = "BDC",
    collection = "SENTINEL-2-16D",
    raster_cube = cube_012014_local,
    vector_dir = segments_path,
    vector_band = "segments",
    version = "snic-30-05",
    parse_info =  c( "X1", "X2","tile", "start_date", "end_date", "band", "version")
)

# Recover random forest model
rf_model <- readRDS("~/PRODES-SITS/models/rf_model_13bands.rds")

# classify vector cube with bands and indices
cube_probs_seg <- sits_classify(
    data = local_segs_cube,
    ml_model = rf_model,
    memsize = 14,
    multicores = 4,
    output_dir = segments_path
)
# obtain class vector cube
cube_class_seg <- sits_label_classification(
    cube_probs_seg,
    memsize = 12,
    multicores = 4,
    output_dir = segments_path
)

# --- Part 6. Visualization

# rename labels to match those with names in the SITS colors
new_labels <- c("Water", "Degraded_Forest",
                "Degradation_Fire", "Clear_Cut_Trees",
                "Clear_Cut_Soil", "Clear_Cut_Vegetation",
                "Forest", "Non_Forest", "Wetlands")

sits_labels(cube_class_seg) <- new_labels

# recover local vector cube
class_seg_local <- sits_cube(
    source = "BDC",
    collection = "SENTINEL-2-16D",
    raster_cube = cube_012014_local,
    vector_dir = segments_path,
    vector_band = "class",
    version = "v1",
    parse_info =  c( "X1", "X2","tile", "start_date", "end_date", "band", "version")
)

# plot the segments cube
plot(class_seg_local)
