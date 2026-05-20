# Load required libraries
library(sf)
library(segmetric)
library(furrr)
library(future)
library(parallelly)
library(openxlsx)

# ============================================================
# Define parameters and paths
# ============================================================

tile <- "012014"
version <- "camara-agrupado-reclassificado"

prodes_dir  <- "data/raw/samples/prodes-2025"
class_dir   <- paste0("data/class/",tile,"/post_processed")

# List the classified polygons 
class <- list.files(
  class_dir,
  pattern    = paste0("^class-post-processed_", tile, "_[^_]+_[^_]+_", version, "\\.gpkg$"),
  full.names = TRUE,
  recursive  = TRUE
)

# List the PRODES reference polygons 
ref_prodes <- list.files(
  prodes_dir,
  pattern    = paste0(".*", tile, ".*\\.gpkg$"),
  full.names = TRUE,
  recursive  = TRUE
)

# ============================================================
# Define intersection metrics function
# ============================================================

compute_metrics <- function(ref_prodes,
                                class,
                                n_workers = floor(parallelly::availableCores() / 2)) {
  
  # --- Configure parallelization ----------------------------
  plan(multisession, workers = n_workers) # Configure the execution plan in parallel using multiple sessions
  on.exit(plan(sequential), add = TRUE)
  message(sprintf("[PARALLEL] Using %d workers.", n_workers))
  
  old_s2_state <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE) # Disables spherical geometry (S2) from the sf package, using planar geometry (GEOS) instead
  on.exit(sf::sf_use_s2(old_s2_state), add = TRUE) # Ensures that spherical geometry is reactivated when this auxiliary function finishes
  
  # --- Helper: secure reading with correction via GEOS --------
  safe_read_sf <- function(path, layer_name) { # Creates an internal function to read spatial files
    layer <- sf::read_sf(path) # Reads the spatial file at the specified path and stores it in the 'layer' object
    
    layer <- sf::st_make_valid(layer) # Corrects topology via GEOS (works on lon/lat with s2 disabled)

    layer <- sf::st_collection_extract(layer, "POLYGON")
    
    valid_mask <- sf::st_is_valid(layer) & !sf::st_is_empty(layer) # Remove geometries that are still invalid or empty
    n_dropped  <- sum(!valid_mask)
    if (n_dropped > 0){
      warning(sprintf("[%s] %d Invalid geometry(ies) removed after correction.",
                      layer_name, n_dropped))
    
    layer <- layer[valid_mask, ]}
    
    return(layer)
  }
 
  # --- Load PRODES (Reference) -------------------------------------
  # Checks if the user passed more than one reference file
  if (length(ref_prodes) > 1)
    warning(sprintf("[PRODES] %d files found; using the first one.", length(ref_prodes)))

  # Read and correct the first PRODES file using the auxiliary function created above
  pol_ref <- safe_read_sf(ref_prodes[1], "PRODES")
  
  # --- Load SITS classification -------------------------
  # Checks if the user has passed more than one classification file
  if (length(class) > 1)
    warning(sprintf("[SITS] %d files found; using the first one.", length(class)))

  # Read and correct the first classification file using the auxiliary function created above
  pol_class <- safe_read_sf(class[1], "SITS")
  
  # --- Align CRS -----------------------------------------
  if (!sf::st_crs(pol_ref) == sf::st_crs(pol_class)) {
    message("[CRS] Reprojecting pol_class for the CRS of pol_ref.")
    pol_class <- sf::st_transform(pol_class, sf::st_crs(pol_ref))
  }
  
  # --- Segmetric Object -----------------------
  # Creates a base object from the 'segmetric' package by crossing the reference polygons and segmentation/classification
  seg_obj <- segmetric::sm_read(ref_sf = pol_ref, seg_sf = pol_class)
  
  # --- Calculate metrics in parallel (segmetric) -----------
  metricas_ids <- c(
    "OS2", "US2", "AFI", "D_index",
    "precision", "recall", "M", "IoU", "Dice"
  )

  # Maps (iterates) the calculation of each metric in parallel using the 'furrr' package
  metricas_list <- furrr::future_map(
    metricas_ids,
    .f        = \(m) {sf::sf_use_s2(FALSE)
                      list(name = m,
                           result = segmetric::sm_compute(seg_obj, m))}, # For each metric 'm', create a list with the metric name and the result of the calculation via segmetric
    .options  = furrr::furrr_options(seed = TRUE), # Ensures reproducibility with reliable parallel generation of random numbers
    .progress = TRUE
  )
  
  # --- Summarize metrics in data.frame (structure: 1 line per metric) ----
  r_df <- do.call(rbind, lapply(metricas_list, function(x) { # Combines the list of results into a single data frame by joining the rows
    
    # Extracts only the numeric values from the segmetric object (remove infinites (Inf) or NA values)
    vals <- unlist(x$result[sapply(x$result, is.numeric)])
    vals <- vals[is.finite(vals)]   # remove NA/Inf

    # Creates a data frame row with descriptive statistics for each metric's values
    data.frame(
      Minimo        = min(vals,            na.rm = TRUE),
      Maximo        = max(vals,            na.rm = TRUE),
      Media         = mean(vals,           na.rm = TRUE),
      Mediana       = median(vals,         na.rm = TRUE),
      Desvio_Padrao = sd(vals,             na.rm = TRUE),
      row.names     = x$name
    )
  }))
  
  # --- Calculate the area of intersection -------------------------
  
  # Total area of PRODES e total area of classification
  area_ref_ha   <- sum(as.numeric(sf::st_area(pol_ref)),   na.rm = TRUE) / 10000
  area_class_ha <- sum(as.numeric(sf::st_area(pol_class)), na.rm = TRUE) / 10000
  
  intersec_sf  <- sf::st_intersection(pol_ref, pol_class)
  areas_m2     <- as.numeric(sf::st_area(intersec_sf))
  areas_ha     <- areas_m2 / 10000
  total_area_ha <- sum(areas_ha, na.rm = TRUE)
  
  # Percentage of intersection relative to PRODES → geometric recall
  pct_ref   <- (total_area_ha / area_ref_ha)   * 100
  
  # Percentage of intersection relative to classification → geometric precision
  pct_class <- (total_area_ha / area_class_ha) * 100
  
  area_df <- data.frame(
    polygon_id      = seq_along(areas_ha),
    area_intersecao_ha = areas_ha
  )
  
  summary_area_df <- data.frame(
    descricao = c(
      "PRODES area (ha)",
      "Classification area (ha)",
      "Intersection area (ha)",
      "Intersection / PRODES (%)",
      "Intersection / Classificacao (%)"
    ),
    valor = c(
      area_ref_ha,
      area_class_ha,
      total_area_ha,
      pct_ref,
      pct_class
    )
  )
  
  message(sprintf("Total area of PRODES:  %.4f ha", area_ref_ha))
  message(sprintf("Total area of classification: %.4f ha", area_class_ha))
  message(sprintf("Total area of intersection: %.4f ha", total_area_ha))
  message(sprintf("intersection / PRODES:  %.2f %%", pct_ref))
  message(sprintf("intersection / Class.:  %.2f %%", pct_class))

  # --- Export to Excel ---------------------------------
  class_filename <- tools::file_path_sans_ext(basename(class[1])) # Extracts the original sorting filename, without the full path and without the extension
  metrics_filename <- sub("^class", "metrics", class_filename) # Replace the "class" prefix (if it exists at the beginning of the string) with "metrics" to generate the new filename
  
  wb <- openxlsx::createWorkbook() # Creates a new "Workbook" (blank Excel file in memory) from the openxlsx package

  # Extract the name of the PRODES file (without extension) to use as the name of the first spreadsheet tab, add a new tab to Excel with the name of the reference file, and write the data frame with the statistics for the 'r_df' metrics in the created tab
  sheet_name <- tools::file_path_sans_ext(basename(ref_prodes[1])) 
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet_name, r_df, rowNames = TRUE)

  # Add a second tab to Excel called "Polygons_Intersection" and write the data frame with the detailed area of ​​each polygon in the second tab
  openxlsx::addWorksheet(wb, "Polygons_Intersection") 
  openxlsx::writeData(wb, "Polygons_Intersection", area_df) # 

  # Add a third tab called "Summary_Intersection" and write the data frame with the summary of total areas and percentages in the third tab
  openxlsx::addWorksheet(wb, "Summary_Intersection") 
  openxlsx::writeData(wb, "Summary_Intersection", summary_area_df)

  # Defines the full path where the Excel file will be saved (in the same directory as the input sorting file)
  out_path <- file.path(
    dirname(class[1]), # Get the folder where the 'class' file is located
    paste0(metrics_filename, ".xlsx") # Attach the Excel file name
  )
  
  # It physically saves the Excel file to the hard drive, overwriting any that already exist
  openxlsx::saveWorkbook(wb, out_path, overwrite = TRUE)
  message(sprintf("[OK] File saved in: %s", out_path))

  # Returns a list in 'invisible' form
  invisible(list(
    metricas        = r_df,
    area_intersecao = area_df,
    total_area_ha      = total_area_ha,
    out_path        = out_path
  ))
}

# ============================================================
# Use intersection metrics function
# ============================================================

metrics_results <- compute_metrics(
  ref_prodes = ref_prodes,
  class   = class
)
