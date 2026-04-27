# ============================================================
# Post-Processing in Forest Areas
# ============================================================

# ============================================================
# 1. Libraries, paths, and initial parameters
# ============================================================

# Step 1.1 -- Load required libraries
library(sf)
library(dplyr)
library(sits)
library(terra)
library(units)
library(smoothr)

# Step 1.2 -- Define paths for files and folders
tile            <- "012014"
version         <- "rf-2y-novos-segmentos"
class_path      <- "data/class"
mask_path       <- "data/raw/auxiliary/masks" #nome da máscara em gpkg geral


post_class_path <- file.path(class_path, tile, "post_processed")
dir.create(post_class_path,
           showWarnings = FALSE,
           recursive = TRUE)

pattern <- paste0(".*", tile, ".*", "class",
                  ".*", version, ".*\\.gpkg$")

raw_class_path <- list.files(class_path,
                             pattern = pattern,
                             full.names = TRUE)

raw_class <- read_sf(raw_class_path)

# ============================================================
# 2. Probabilistic reclassification
# ============================================================

# Step 2.1 -- Create sum columns 
raw_class <- raw_class |>
  mutate(
    Suppression_sum = rowSums(
      across(
        all_of(c(
          "Corte_Raso",
          "Corte_Raso_Com_Vegetacao",
          "Corte_Raso_Com_Fogo",
          "Corte_Raso_Com_Arvores_Remanescentes"
        ))
      ),
      na.rm = TRUE
    ),
    Degrad_sum = rowSums(
      across(
        all_of(c(
          "Degradacao",
          "Degradacao_Por_Fogo"
        ))
      ),
      na.rm = TRUE
    ),
    Natural = rowSums(
      across(
        all_of(c(
          "Floresta",
          "Floresta_Transicional",
          "Vegetacao_Natural_Nao_Florestal"
        ))
      ),
      na.rm = TRUE
    )
  )

# Step 2.4 -- Reclassification
raw_class <- raw_class |>
  mutate(
    class = if_else(
      Suppression_sum >= Degrad_sum &
        Suppression_sum >= Natural &
        Suppression_sum >= Corpo_Dagua,
      "Soma_Supressao",
      class
    )
  )

# Step 2.5 -- Filter only suppression-related classes
post_class <- raw_class |>
  filter(
    class %in% c(
      "Corte_Raso",
      "Corte_Raso_Com_Vegetacao",
      "Corte_Raso_Com_Fogo",
      "Corte_Raso_Com_Arvores_Remanescentes",
      "Suppression_sum"
    )
  ) |>
  mutate(class = "supressao")

# Step 2.7 -- Dissolve and aggregate geometries
post_class <- post_class |>
  group_by(class) |>
  summarise(
    geometry = st_union(geom),
    .groups = "drop"
  ) |>
  st_cast("MULTIPOLYGON")

# ============================================================
# 3. Extraction of cloud features
# ============================================================

# Step 3.1 -- Function definition
extract_cloud_mask <- function(
    sits_classification_path,
    sits_reclassification,
    cloud_values = c(3, 8, 9, 10),
    output_dir = NULL,
    collection = "SENTINEL-2-16D",
    date_window_days = 1
) {
  
  # ----------------------------------------------------------
  # 1. Extract metadata from filename
  # ----------------------------------------------------------
  filename_base <- basename(sits_classification_path)
  
  last_date_str <- regmatches(
    filename_base,
    gregexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", filename_base)
  )[[1]]
  
  if (length(last_date_str) == 0) {
    stop("No date in YYYY-MM-DD format found in file name: ", filename_base)
  }
  
  last_date <- as.Date(tail(last_date_str, 1))
  
  start_date_scl <- format(last_date - date_window_days, "%Y-%m-%d")
  end_date_scl   <- format(last_date, "%Y-%m-%d")
  
  tile_id <- regmatches(
    filename_base,
    regexpr("(?<=SENTINEL-2_MSI_)[0-9]+", filename_base, perl = TRUE)
  )
  
  if (length(tile_id) == 0 || tile_id == "") {
    stop("Could not extract tile_id from file name: ", filename_base)
  }
  
  message(" -> Tile: ", tile_id)
  message(" -> SCL time window: ", start_date_scl, " to ", end_date_scl)
  
  # ----------------------------------------------------------
  # 2. Build BDC cube
  # ----------------------------------------------------------
  scl_cube <- sits::sits_cube(
    source = "BDC",
    collection = collection,
    tiles = tile_id,
    bands = "SCL",
    start_date = start_date_scl,
    end_date = end_date_scl
  )
  
  # ----------------------------------------------------------
  # 3. Load SCL raster
  # ----------------------------------------------------------
  scl_files <- scl_cube$file_info[[1]] |>
    dplyr::filter(
      band == "CLOUD",
      date == as.Date(end_date_scl)
    ) |>
    dplyr::pull(path)
  
  if (length(scl_files) == 0) {
    stop("No SCL file found for date: ", end_date_scl)
  }
  
  scl_raster <- terra::rast(scl_files[1])
  
  message(" -> SCL file loaded: ", scl_files[1])
  
  # ----------------------------------------------------------
  # 4. Create cloud mask
  # ----------------------------------------------------------
  scl_mask <- terra::classify(
    scl_raster,
    rcl = cbind(cloud_values, rep(1, length(cloud_values))),
    others = NA
  )
  
  # ----------------------------------------------------------
  # 5. Crop to classification extent
  # ----------------------------------------------------------
  class_bbox <- sits_reclassification |>
    sf::st_transform(terra::crs(scl_mask)) |>
    terra::vect() |>
    terra::ext()
  
  scl_raster_crop <- terra::crop(scl_mask, class_bbox)
  
  # ----------------------------------------------------------
  # 6. Vectorize cloud mask
  # ----------------------------------------------------------
  cloud_vec <- terra::as.polygons(scl_raster_crop, dissolve = TRUE) |>
    sf::st_as_sf() |>
    sf::st_transform(sf::st_crs(sits_reclassification)) |>
    smoothr::fill_holes(threshold = Inf)
  
  # ----------------------------------------------------------
  # 7. Save (optional)
  # ----------------------------------------------------------
  if (!is.null(output_dir)) {
    output_filename <- paste0(
      "cloud-vec_", tile_id, "_", end_date_scl, ".gpkg"
    )
    
    output_path <- file.path(output_dir, output_filename)
    
    cloud_vec |>
      sf::st_transform(4674) |>
      sf::st_write(output_path, append = FALSE)
    
    message(" -> Cloud vector saved (EPSG:4674): ", output_path)
  }
  
  # ----------------------------------------------------------
  # 8. Return
  # ----------------------------------------------------------
  return(
    invisible(
      list(
        cloud_vec   = cloud_vec,
        tile_id     = tile_id,
        end_date_scl = end_date_scl
      )
    )
  )
}

# Step 3.2 -- Apply function
result <- extract_cloud_mask(
  sits_classification_path = raw_class_path,
  sits_reclassification    = post_class,
  cloud_values             = c(3, 8, 9, 10),
  output_dir               = post_class_path
)

# Step 3.3 -- Extract outputs
cloud_vec   <- result$cloud_vec
end_date_scl <- result$end_date_scl

# ============================================================
# 4. Cloud/shadow difference
# ============================================================

Cloud_union <- sf::st_union(cloud_vec)

cloud_vec_buffer <- sf::st_buffer(Cloud_union, dist = 100)

reclass_cloud_cleaned <- sf::st_difference(
        post_class,
        cloud_vec_buffer
      ) |>
      sf::st_cast("MULTIPOLYGON")

# ============================================================
# 5. Fill holes < 1 hectare (first round)
# ============================================================

query <- sprintf("SELECT * FROM  WHERE tile = '%s'", tile) #nome da máscara em gpkg geral
prodes_mask <- read_sf(mask_path,
                       query = query) 

prodes_mask <- sf::st_transform(
  prodes_mask,
  sf::st_crs(reclass_cloud_cleaned)
)

merged <- list(reclass_cloud_cleaned, prodes_mask) |>
  purrr::map(sf::st_make_valid) |>
  purrr::map(\(x) sf::st_transform(x, sf::st_crs(reclass_cloud_cleaned))) |>
  purrr::map(\(x) {
    sf::st_geometry(x) <- "geom"
    x
  }) |>
  purrr::map(\(x) sf::st_cast(x, "MULTIPOLYGON")) |>
  dplyr::bind_rows() |>
  sf::st_union()

smoothed <- smoothr::fill_holes(
  merged,
  threshold = units::set_units(10000, "m^2")
)

# ============================================================
# 6. Difference with deforestation mask (first round)
# ============================================================

smoothed <- sf::st_transform(smoothed, sf::st_crs(prodes_mask))

smoothed <- smoothed |> 
  st_make_valid()

prodes_mask <- prodes_mask |>
  st_make_valid()

mask_union <- prodes_mask |>
  st_union() |>
  st_make_valid()

class_diff_mask <- sf::st_difference(
  smoothed,
  mask_union
) |>
  sf::st_cast("POLYGON") |>
  sf::st_sf()

# ============================================================
# 7. Fill bays and smooth edges 
# ============================================================

class_diff_mask_filled_bays <- class_diff_mask |>
  sf::st_buffer(dist = 50) |>
  sf::st_buffer(dist = -50) |>
  sf::st_cast("POLYGON") |>
  tibble::rowid_to_column("id")

sf::st_write(
  class_diff_mask_filled_bays,
  file.path(
    output_dir,
    paste0("class_diff_mask_filled_bays_", tile_id, "_", end_date_scl, ".gpkg")
  )
)
# vai salvar???????
# ============================================================
# 8. Fill holes < 1 hectare (second round)
# ============================================================

merged_2 <- list(class_diff_mask_filled_bays, mask) |>
  purrr::map(sf::st_make_valid) |>
  purrr::map(\(x) sf::st_transform(x, sf::st_crs(class_diff_mask_filled_bays))) |>
  purrr::map(\(x) {
    sf::st_geometry(x) <- "geom"
    x
  }) |>
  purrr::map(\(x) sf::st_cast(x, "MULTIPOLYGON")) |>
  dplyr::bind_rows() |>
  sf::st_union()

smoothed_2 <- smoothr::fill_holes(
  merged_2,
  threshold = units::set_units(10000, "m^2")
)

sf::st_write(
  smoothed_2,
  file.path(
    output_dir,
    paste0("smoothed_2_", tile_id, "_", end_date_scl, ".gpkg")
  )
)
# vai salvar???????
# ============================================================
# 9. Difference with deforestation mask (second round)
# ============================================================

class_diff_mask_2 <- sf::st_difference(
  smoothed_2,
  mask_union
) |>
  sf::st_cast("POLYGON") |>
  sf::st_sf()

sf::st_write(
  class_diff_mask_2,
  file.path(
    output_dir,
    paste0("class_diff_mask_2_", tile_id, "_", end_date_scl, ".gpkg")
  )
)
#pra que recotar pela máscara de novo?????

# ============================================================
# 10. Remove polygons < 1 hectare
# ============================================================

class_diff_mask_2$area_m2 <- as.numeric(sf::st_area(class_diff_mask_2))
class_diff_mask_2$area_ha <- class_diff_mask_2$area_m2 / 10000

class_diff_mask_bigger_than_1ha <- class_diff_mask_2 |>
  dplyr::filter(area_ha >= 1)

# ============================================================
# 11. Save final result
# ============================================================

poligonos_supressao <- st_transform(
  class_diff_mask_bigger_than_1ha,
  crs = 4674
) |>
  sf::st_cast("POLYGON") 

sf::st_write(
  poligonos_supressao,
  file.path(
    post_class_path,
    paste0("sits-classification-post-processed_",
           tile, "_", end_date_scl, ".gpkg")
  )
)

# ============================================================
# 12. Spatial logic (adapted for sf)
# ============================================================

# Step A: Separate classes
arvore_remanesce <- result_filtered %>%
  filter(class == "DESMAT_ARVORE_REMANESCE")

outras_classes <- result_filtered %>%
  filter(class != "DESMAT_ARVORE_REMANESCE")

if (nrow(arvore_remanesce) > 0 && nrow(outras_classes) > 0) {
  
  message("Calculating spatial relationships...")
  
  # Step B: Reference geometry
  outras_classes_union <- st_union(outras_classes)
  
  # Step C: Disjoint test
  matriz_disjoint <- st_disjoint(
    arvore_remanesce,
    outras_classes_union,
    sparse = FALSE
  )
  
  is_separated <- matriz_disjoint[, 1]
  
  # Step D: Filter (keep intersecting)
  arvore_remanesce_limpas <- arvore_remanesce[!is_separated, ]
  
  # Step E: Merge results
  final_result <- bind_rows(
    outras_classes,
    arvore_remanesce_limpas
  )
  
} else {
  final_result <- result_filtered
}

