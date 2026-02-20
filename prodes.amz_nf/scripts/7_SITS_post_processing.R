library(sf)
library(dplyr)

# ============================================================
# INPUT PARAMETERS
# ============================================================

class_obj           <- "SENTINEL-2_MSI_012014_2024-07-27_2025-07-28_class_1y-012014"
classification_path <- paste0("data/class/", class_obj, ".gpkg") 
prodes_mask_path    <- "data/raw/auxiliary/grade_bdc_acumulado_24_10857.gpkg"
output_path         <- paste0("data/class/post_process/", class_obj, "_post-processed.gpkg")

tile <- stringr::str_split_i(class_obj, "-", -1) 

crs_albers <- "+proj=aea +lat_0=-12 +lon_0=-54 +lat_1=-2 +lat_2=-22 +x_0=5000000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"

# ============================================================
# PROCESSING
# ============================================================

# 1. Load data
classification <- read_sf(classification_path) 
prodes_mask    <- read_sf(prodes_mask_path) 

# 2. Filter deforestation classes
deforestation_classes <- c(
  "DESMAT_ARVORE_REMANESCE",
  "DESMAT_CORTE_RASO",
  "DESMAT_DEGRAD_FOGO",
  "DESMAT_VEG"
)

classification_filtered <- classification %>% 
  filter(class %in% deforestation_classes) %>% 
  st_transform(crs_albers)

# 3. Dissolver filtered polygons & Desagregar (Explodir multipartes)
classification_dissolved <- classification_filtered %>%
  group_by(class) %>%
  summarise(geom = st_union(geom)) %>% # Agrega (Dissolve)
  st_cast("MULTIPOLYGON") %>%          # Garante tipo uniforme
  st_cast("POLYGON")                   # Desagrega (disagg): Transforma Multipart em Singlepart

# 4. Extract PRODES mask by tile
mask_tile <- prodes_mask %>% 
  filter(tile == !!tile) %>% 
  st_transform(crs_albers)

# 5. Fix mask geometries
mask_tile <- st_make_valid(mask_tile)

# 6. Difference (classification - mask)
message("Aplicando máscara (difference)...")

# st_difference remove a intersecção
difference <- st_difference(classification_dissolved, mask_tile)

difference <- difference %>%
  st_make_valid() %>%
  st_collection_extract("POLYGON") %>% # Garante que não sobrou linha/ponto
  st_cast("POLYGON") # Garante polígonos simples (disagg)

# 9. Calculate area in hectares.
difference <- difference %>%
  mutate(area_ha = as.numeric(st_area(geom)) / 10000) %>% 
  ungroup() %>%                      # 1. Garante que não há grupos presos
  filter(!st_is_empty(.)) %>%        # 2. REMOVE GEOMETRIAS VAZIAS
  st_make_valid()

# 10. Extract features with area >= 1 ha
result_filtered <- difference %>%
  filter(area_ha >= 1)

# ============================================================
# 13. LÓGICA ESPACIAL (ADAPTADA PARA SF)
# ============================================================

# Passo A: Separar
arvore_remanesce <- result_filtered %>% filter(class == "DESMAT_ARVORE_REMANESCE")
outras_classes   <- result_filtered %>% filter(class != "DESMAT_ARVORE_REMANESCE")

if (nrow(arvore_remanesce) > 0 && nrow(outras_classes) > 0) {
  
  message("Calculando relações espaciais...")
  
  # Passo B: Criar a geometria de referência (União das outras classes)
  outras_classes_union <- st_union(outras_classes)
  
  # Passo C: Verificar disjoint (Separado)
  # sparse = FALSE retorna uma matriz densa (TRUE/FALSE) em vez de lista de índices
  # st_disjoint retorna TRUE se ESTIVER SEPARADO
  matriz_disjoint <- st_disjoint(arvore_remanesce, outras_classes_union, sparse = FALSE)
  
  # Transformar em vetor lógico (coluna 1, pois união é 1 geometria só)
  is_separated <- matriz_disjoint[, 1]
  
  # Passo D: Filtrar
  # Queremos MANTER as árvores que NÃO estão separadas (ou seja, tocam/interceptam)
  # Lógica: Se is_separated é FALSE (toca), !is_separated é TRUE (mantém)
  arvore_remanesce_limpas <- arvore_remanesce[!is_separated, ]
  
  # Passo E: Juntar de volta (rbind do sf é igual ao do base R, ou bind_rows do dplyr)
  final_result <- bind_rows(outras_classes, arvore_remanesce_limpas)
  
} else {
  final_result <- result_filtered
}

# ============================================================
# OUTPUT
# ============================================================

arvore_remanesce <- arvore_remanesce %>%
  filter(area_ha >= 6.25)
final_result <- bind_rows(outras_classes, arvore_remanesce)
                          
st_write(final_result, output_path, delete_dsn = TRUE)