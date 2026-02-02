# ============================================================
# Classification Post-processing Script with PRODES Mask
# ============================================================

library(sf)
library(dplyr)

# ============================================================
# INPUT PARAMETERS
# ============================================================

# Define your data paths
classification_path <- "path/to/classification.shp"
prodes_mask_path <- "path/to/prodes_mask.shp"
output_path <- "path/to/output/classification_post_processed.shp"

# Parameters
tile <- "012014"  # Tile to process
dissolve_by_class <- TRUE  # Dissolve by deforestation class?

# ============================================================
# PROCESSING
# ============================================================

# 1. Load data
classification <- st_read(classification_path, quiet = TRUE)
prodes_mask <- st_read(prodes_mask_path, quiet = TRUE)

# 2. Extract by expression - Filter deforestation classes
deforestation_classes <- c(
  "DESMAT_ARVORE_REMANESCE",
  "DESMAT_CORTE_RASO",
  "DESMAT_DEGRAD_FOGO",
  "DESMAT_VEG",
  "DESMAT_CORTE_RASO_DM",
  "DESMAT_VEG_DM"
)

classification_filtered <- classification %>%
  filter(class %in% deforestation_classes)

# 3. Dissolve geometries
if (dissolve_by_class) {
  # Dissolve by class
  classification_dissolved <- classification_filtered %>%
    group_by(class) %>%
    summarise(geometry = st_union(geometry), .groups = "drop")
} else {
  # Dissolve everything into a single geometry
  classification_dissolved <- classification_filtered %>%
    summarise(geometry = st_union(geometry))
}

# 4. Extract PRODES mask by tile
mask_tile <- prodes_mask %>%
  filter(tile == !!tile)

# 5. Fix mask geometries
mask_tile <- st_make_valid(mask_tile)

# 6. Difference (classification - mask)
difference <- st_difference(classification_dissolved, st_union(mask_tile))

# 7. Convert multiparts to single parts
difference_single <- st_cast(difference, "POLYGON")

# 8. Drop 'fid' field if exists
if ("fid" %in% names(difference_single)) {
  difference_single <- difference_single %>%
    select(-fid)
}

# 9. Calculate area in hectares
difference_single <- difference_single %>%
  mutate(area_ha = as.numeric(st_area(geometry)) / 10000)

# 10. Extract features with area >= 1 ha
result_filtered <- difference_single %>%
  filter(area_ha >= 1)

# 11. Promote to multiparts
result_multi <- result_filtered %>%
  group_by(across(-geometry)) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_cast("MULTIPOLYGON")

# 12. Extract only DESMAT_ARVORE_REMANESCE
if ("class" %in% names(result_multi)) {
  desmat_tree <- result_multi %>%
    filter(class != "DESMAT_ARVORE_REMANESCE")
} else {
  desmat_tree <- result_multi
}

# 13. Extract by location (features that DO NOT intersect or touch DESMAT_ARVORE_REMANESCE)
if (nrow(desmat_tree) > 0 && "class" %in% names(result_multi)) {
  # Create union of DESMAT_ARVORE_REMANESCE geometries
  desmat_tree_union <- result_multi %>%
    filter(class == "DESMAT_ARVORE_REMANESCE") %>%
    st_union()
  
  # Check which features DO NOT intersect
  if (length(desmat_tree_union) > 0) {
    intersects <- st_intersects(result_multi, desmat_tree_union, sparse = FALSE)
    touches <- st_touches(result_multi, desmat_tree_union, sparse = FALSE)
    
    # Keep only features that DO NOT intersect AND DO NOT touch
    final_result <- result_multi[!intersects[,1] & !touches[,1], ]
  } else {
    final_result <- result_multi
  }
} else {
  final_result <- result_multi
}

# 14. Save result
st_write(final_result, output_path, delete_dsn = TRUE, quiet = TRUE)