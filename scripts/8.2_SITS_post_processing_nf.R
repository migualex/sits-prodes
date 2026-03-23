# =============================================================================
# STAGE 0 — Readings and initial preparation
# =============================================================================

library(sf)
library(dplyr)
library(sits)
library(terra)
library(units)
library(smoothr)

# Read SITS classification
mask_path <- "~/grupos/biomasbr/amazonia/sits-prodes/prodes.amz_nf/data/raw/auxiliary/mask_geral_nf.gpkg"
sits_classification_path <- "~/grupos/biomasbr/amazonia/sits-prodes/prodes.amz_nf/data/class/015002/2y/SENTINEL-2_MSI_015002_2023-08-13_2025-07-28_class_rf-2y-015002-nf-classes-novas.gpkg"
sits_classification_raw <- sf::st_read(sits_classification_path, quiet = TRUE)
output_dir <- sub("/SENTINEL.*", "", sits_classification_path)

# =============================================================================
# STAGE 1 — Probabilistic reclassification
# =============================================================================

# Step 1.1 -- Create column "Soma_Supressao"
sits_classification_raw <- sits_classification_raw |>
  mutate(
    Soma_Supressao = rowSums(across(all_of(c(
      "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Agricultura",
      "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Solo_Exposto",
      "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Agricultura_Antigo",
      "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Solo_Exposto_Antigo"
    ))), na.rm = TRUE)
  )

# Step 1.2 -- Create column "Soma_NF"
sits_classification_raw <- sits_classification_raw |>
  mutate(
    Soma_NF = rowSums(across(all_of(c(
      "Vegetacao_Natural_Nao_Florestal_Herbacea_Seca_Mais_Biomassa",
      "Vegetacao_Natural_Nao_Florestal_Herbacea_Seca_Menos_Biomassa",
      "Vegetacao_Natural_Nao_Florestal_Herbacea_Umida",
      "Vegetacao_Natural_Nao_Florestal_Vereda",
      "Vegetacao_Natural_Nao_Florestal_Mata",
      "Vegetacao_Natural_Nao_Florestal_Transicao_Florestal",
      "Vegetacao_Natural_Nao_Florestal_Herbacea_Seca_Pos_Fogo"
    )
    )), na.rm = TRUE)
  )

# Step 1.3 -- Create column "Soma_Agua"
sits_classification_raw <- sits_classification_raw |>
  mutate(
    Soma_Agua = rowSums(across(all_of(c(
      "Hidrografia_Lago",
      "Hidrografia_Rio"
    ))), na.rm = TRUE)
  )

# Step 1.4 -- Reclassify all segments where Sum_Suppression > Sum_NF & Sum_Suppression >= Sum_Water
sits_classification_raw <- sits_classification_raw |>
  mutate(
    class = if_else(Soma_Supressao >= Soma_NF & Soma_Supressao >= Soma_Agua, "Soma_Supressao", class)
  )

# Step 1.5 -- Filter the dataset to keep only the classes related to suppression
sits_classification <- sits_classification_raw |>
  dplyr::filter(class %in% c(
    "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Agricultura",
    "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Solo_Exposto",
    "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Agricultura_Antigo",
    "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Solo_Exposto_Antigo",
    "Soma_Supressao" # adicionado no passo anterior
  ))

# Step 1.6 -- Standardize the nomenclature of the classes
sits_reclassification <- sits_classification |>
  dplyr::mutate(class = "supressao")

# Step 1.7 -- Dissolving and aggregating geometries
sits_reclassification <- sits_reclassification |>
  dplyr::group_by(class) |>
  dplyr::summarise(geometry = sf::st_union(geom), .groups = "drop") |>
  sf::st_cast("MULTIPOLYGON")

# Salvar resultado parcial
sf::st_write(sits_reclassification, file.path(output_dir, "sits_reclassification.gpkg"))

# =============================================================================
# STAGE 2 — Extraction of cloud features
# 
# Extracts vector cloud/shadow mask from the SCL band at the latest date
# from the BDC cube that generated the sits classification
# =============================================================================

# Step 2.1 -- Defining the 'extract_cloud_mask' function
extract_cloud_mask <- function(
    sits_classification_path,
    sits_reclassification,
    cloud_values    = c(3, 8, 9, 10),
    output_dir      = NULL,
    collection      = "SENTINEL-2-16D",
    date_window_days = 1
) {
  # ---------------------------------------------------------------------------
  # 1. Extrair metadados do nome do arquivo
  # ---------------------------------------------------------------------------
  filename_base <- basename(sits_classification_path)
  
  last_date_str <- regmatches(
    filename_base,
    gregexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", filename_base)
  )[[1]]
  
  if (length(last_date_str) == 0) {
    stop("Nenhuma data no formato YYYY-MM-DD encontrada no nome do arquivo: ", filename_base)
  }
  
  last_date      <- as.Date(tail(last_date_str, 1))
  start_date_scl <- format(last_date - date_window_days, "%Y-%m-%d")
  end_date_scl   <- format(last_date, "%Y-%m-%d")
  
  tile_id <- regmatches(
    filename_base,
    regexpr("(?<=SENTINEL-2_MSI_)[0-9]+", filename_base, perl = TRUE)
  )
  
  if (length(tile_id) == 0 || tile_id == "") {
    stop("Não foi possível extrair o tile_id do nome do arquivo: ", filename_base)
  }
  
  message("  -> Tile: ", tile_id)
  message("  -> Janela temporal SCL: ", start_date_scl, " a ", end_date_scl)
  
  # ---------------------------------------------------------------------------
  # 2. Montar cubo BDC com banda SCL
  # ---------------------------------------------------------------------------
  scl_cube <- sits::sits_cube(
    source     = "BDC",
    collection = collection,
    tiles      = tile_id,
    bands      = "SCL",
    start_date = start_date_scl,
    end_date   = end_date_scl
  )
  
  # ---------------------------------------------------------------------------
  # 3. Carregar raster SCL
  # ---------------------------------------------------------------------------
  scl_files <- scl_cube$file_info[[1]] |>
    dplyr::filter(band == "CLOUD", date == as.Date(end_date_scl)) |>
    dplyr::pull(path)
  
  if (length(scl_files) == 0) {
    stop("Nenhum arquivo SCL encontrado para a data: ", end_date_scl)
  }
  
  scl_raster <- terra::rast(scl_files[1])
  message("  -> Arquivo SCL carregado: ", scl_files[1])
  
  # ---------------------------------------------------------------------------
  # 4. Criar máscara binária de nuvens/sombras
  # ---------------------------------------------------------------------------
  scl_mask <- terra::classify(
    scl_raster,
    rcl    = cbind(cloud_values, rep(1, length(cloud_values))),
    others = NA
  )
  
  # ---------------------------------------------------------------------------
  # 5. Recortar pela extensão da classificação
  # ---------------------------------------------------------------------------
  class_bbox <- sits_reclassification |>
    sf::st_transform(terra::crs(scl_mask)) |>
    terra::vect() |>
    terra::ext()
  
  scl_raster_crop <- terra::crop(scl_mask, class_bbox)
  
  # ---------------------------------------------------------------------------
  # 6. Vetorizar máscara de nuvens
  # ---------------------------------------------------------------------------
  cloud_vec <- terra::as.polygons(scl_raster_crop, dissolve = TRUE) |>
    sf::st_as_sf() |>
    sf::st_transform(sf::st_crs(sits_reclassification))
  
  # ---------------------------------------------------------------------------
  # 7. Salvar resultado (opcional)
  # ---------------------------------------------------------------------------
  if (!is.null(output_dir)) {
    output_filename <- paste0(
      "cloud_vec_",
      tile_id, "_",
      end_date_scl,       # formato já é YYYY-MM-DD
      ".gpkg"
    )
    output_path <- file.path(output_dir, output_filename)
    sf::st_write(cloud_vec, output_path, append = FALSE)
    message("  -> Vetor de nuvens salvo em: ", output_path)
  }

  # ---------------------------------------------------------------------------
  # 8. Retornar lista com vetor de nuvens e metadados extraídos
  # --------------------------------------------------------------------------- 
  return(invisible(list(
    cloud_vec    = cloud_vec,
    tile_id      = tile_id,
    end_date_scl = end_date_scl
  )))
}

# Step 2.2 -- Applying the 'extract_cloud_mask' function
result <- extract_cloud_mask(
  sits_classification_path = sits_classification_path,
  sits_reclassification      = sits_reclassification,
  cloud_values             = c(3, 8, 9, 10),
  output_dir               = output_dir
)

# Step 2.3 -- Access each object in the result individually
cloud_vec    <- result$cloud_vec
tile_id      <- result$tile_id
end_date_scl <- result$end_date_scl

# Recarregar vetor salvo
cloud_vec <- sf::st_read(file.path(
  output_dir,
  paste0("cloud_vec_", tile_id, "_", end_date_scl, ".gpkg")
))

# =============================================================================
# STAGE 3 — Difference with cloud/shadow
# =============================================================================

# Remove all regions covered by cloud/shadow from the classified areas (suppress)
sits_classification_cloud_cleaned <- sf::st_difference(
  sits_reclassification,
  sf::st_union(cloud_vec)
) |>
  sf::st_cast("MULTIPOLYGON")

# Salvar resultado parcial
sf::st_write(sits_classification_cloud_cleaned, 
             file.path(
               output_dir,
               paste0("sits_classification_cloud_cleaned_", tile_id, "_", end_date_scl, ".gpkg")
             ))

# =============================================================================
# STAGE 4 — Fill holes < 1.5 hectares
# =============================================================================

# Step 4.1 -- Read accumulated deforestation data, filtering only for the tile under analysis
mask <- sf::st_read(mask_path, quiet = TRUE) |>
  dplyr::filter(tile == tile_id)

# Step 4.2 -- Ensure that the mask is in the same spatial reference system
mask <- sf::st_transform(
  mask,
  sf::st_crs(sits_classification_cloud_cleaned)
)

# Step 4.3 -- Merge reclassification with the accumulated mask, 
# dissolving all geometries into a single continuous set
merged <- sf::st_union(
  dplyr::bind_rows(
    sits_classification_cloud_cleaned,
    mask
  )
)

# Step 4.4 -- Fills internal holes smaller than 1.5 hectares (15000 m²)
smoothed <- smoothr::fill_holes(
  merged, 
  threshold = units::set_units(15000, "m^2")
)


# Salvar resultado parcial
sf::st_write(smoothed, 
             file.path(
               output_dir,
               paste0("smoothed_", tile_id, "_", end_date_scl, ".gpkg")
             ))


# =============================================================================
# STAGE 5 - Difference with the deforestation mask
# =============================================================================

# Step 5.1 -- Ensure consistent CRS prior to space operation
smoothed <- sf::st_transform(smoothed, sf::st_crs(mask))

# Step 5.2 -- Fix geometries
smoothed <- smoothed |>
  st_make_valid()
mask <- mask |>
  st_make_valid()

# Step 5.3 -- Make a union
mask_union <- mask |>
  st_union() |>
  st_make_valid()

# Step 5.4 -- Remove areas already present in the deforestation mask
class_diff_mask <- sf::st_difference(
  smoothed,
  mask_union
) |>
  sf::st_cast("POLYGON") |>
  sf::st_sf()


# Salvar resultado parcial
sf::st_write(class_diff_mask, 
             file.path(
               output_dir,
               paste0("class_diff_mask_", tile_id, "_", end_date_scl, ".gpkg")
             ))

# =============================================================================
# STAGE 6 — Fill in "bays" and smooth edges
# =============================================================================

# Step 6.1 -- Mathematical morphology / Positive buffer then negative buffer
class_diff_mask_filled_bays <- class_diff_mask |>
  sf::st_buffer(dist =  40)                    |>
  sf::st_buffer(dist = -40)                    |>
  sf::st_cast("POLYGON")                       |>
  tibble::rowid_to_column("id")


# Salvar resultado parcial
sf::st_write(class_diff_mask_filled_bays, 
             file.path(
               output_dir,
               paste0("40m_class_diff_mask_filled_bays_", tile_id, "_", end_date_scl, ".gpkg")
             ))

# =============================================================================
# STAGE 7 — Exclude polygons < 1 hectares
# =============================================================================

# Step 7.1 -- Calculate the area of each polygon (in m²)
class_diff_mask_filled_bays$area_m2 <- as.numeric(sf::st_area(class_diff_mask_filled_bays))

# Step 7.2 -- Convert square meters to hectares
class_diff_mask_filled_bays$area_ha <- class_diff_mask_filled_bays$area_m2 / 10000

# Step 7.3 -- Keep only polygons with an area greater than or equal to 1 hectare
class_diff_mask_bigger_than_1ha <- class_diff_mask_filled_bays |>
  dplyr::filter(area_ha >= 1)

# =============================================================================
# ETAPA 8 —  Save final result
# =============================================================================

# Step 8.1 -- Reproject to EPSG:4674 (SIRGAS 2000)
poligonos_supressao <- st_transform(class_diff_mask_bigger_than_1ha, crs = 4674)

# Step 8.2 -- Save final result
sf::st_write(poligonos_supressao, 
             file.path(
               output_dir,
               paste0("sits_classification_final_", tile_id, "_", end_date_scl, ".gpkg")
             ))
