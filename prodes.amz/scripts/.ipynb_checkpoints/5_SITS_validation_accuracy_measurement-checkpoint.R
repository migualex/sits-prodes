## I. Load Required Libraries
library(tibble)
library(sits)
library(sitsdata)
library(terra)

## II. Set tempdir 
tempdir_r <- "~/Gustavo/SITS/sits-prodes/prodes.amz/data"
dir.create(tempdir_r, showWarnings = FALSE)

# ============================================================
# 1. Load probabilities cube
# ============================================================

# Step 1.1 -- Define the classes of the probability cube
labels <- c("1" = "Clear_Cut_Bare_Soil",
            "2" = "Clear_Cut_Burned_Area", 
            "3" = "Mountainside_Forest", 
            "4" = "Forest", 
            "5" = "Riparian_Forest", 
            "6" = "Clear_Cut_Vegetation", 
            "7" = "Water",
            "8" = "Seasonally_Flooded",
            "9" = "Wetland")

# Step 1.2 -- Directory where the data is stored 
data_dir <- system.file("extdata/Rondonia-Class-2022-Mosaic/", package = "sitsdata")

# Step 1.3 -- Create a probability data cube from a file 
probs_datacube_class <- sits_cube(
  source = "MPC",
  collection = "SENTINEL-2-L2A",
  data_dir = data_dir,
  bands = "class",
  labels = labels,
  version = "mosaic"
)

# Step 1.3.1 --  Plot the classification map
plot(probs_datacube_class)


# ============================================================
# 2. Uncertainty
# ============================================================

# Step 2.1 -- Calculate uncertainty vector cube
uncertainty <- sits_uncertainty(
  vector_cube,
  type = "entropy",
  multicores = 8L,
  memsize = 80L,
  output_dir = class_path,
  version = version
  progress = TRUE
  )

# Step 2.2 -- armazena o caminho em file_info
uncertainty_path <- uncertainty$file_info[[0]]$path

# Step 2.3 -- Carrega os polígonos de segmentos com a entropia
uncertainty_polygons <- st_read(uncertainty_path)


# Step 2.3 -- Pegue o caminho do arquivo gerado no passo anterior
file_path <- sits_values(raster_uncert_cube)$file_info[[1]]$path

# Step 2.3 -- Carregue o raster
r <- rast(file_path)


# Step 2.3 -- Rasterizar o campo de entropia
raster_entropy <- rasterize(
  x = vect(poligonos), 
  y = temp_raster, 
  field = "entropy" # "entropy" é o nome padrão da coluna criada pelo sits_uncertainty
  )

# Step 2.3 -- Aplique a escala e converta para Uint16 (INT2U no terra)
# Multiplicamos por 10.000 para manter a precisão
raster_uint16 <- round(raster_entropy * 10000)

# Step 2.3 -- Plot the resulting uncertainty cube
plot(raster_uint16)

# Step 2.3 -- Salve o arquivo final com o tipo de dado correto
writeRaster(
  raster_uint16, 
  filename = file.path(class_path, "entropy_uint16.tif"),
  datatype = "INT2U",  # Este é o código para Uint16
  overwrite = TRUE,
  gdal = c("COMPRESS=LZW") # Compressão adicional para diminuir o tamanho
)


# ============================================================
# 3. Cross-validation of training data
# ============================================================

# Step 3.1 -- Using k-fold validation
rfor_validate_mt <- sits_kfold_validate(
  samples = samples_matogrosso_mod13q1,
  folds = 5,
  ml_method = sits_rfor(),
  multicores = 5
)
# Step 3.2 -- Plot the confusion matrix
plot(rfor_validate_mt, type = "confusion_matrix")

# Step 3.3 -- Plot the metrics by class
plot(rfor_validate_mt, type = "metrics")


# ============================================================
# 4. Accuracy assessment of classified images
# ============================================================

# Step 4.1 -- Sampling design
sampling_design <- sits_sampling_design(
  cube = probs_datacube_class,
  expected_ua = c(
    "Clear_Cut_Bare_Soil" = 0.75,
    "Clear_Cut_Burned_Area" = 0.70, 
    "Mountainside_Forest" = 0.70, 
    "Forest" = 0.75,  
    "Riparian_Forest" = 0.70, 
    "Clear_Cut_Vegetation" = 0.70,  
    "Water" = 0.70, 
    "Seasonally_Flooded" = 0.70, 
    "Wetland" = 0.70
  ),
  alloc_options = c(120, 100),
  std_err = 0.01,
  rare_class_prop = 0.1
)

# Step 4.2 -- Show sampling design
sampling_design

# Step 4.3 -- Generate stratified random samples
samples_sf <- sits_stratified_sampling(
  cube = probs_datacube_class,
  sampling_design = sampling_design,
  alloc = "alloc_120",
  multicores = 4
)

# Step 4.4 -- Save samples in a shapefile
sf::st_write(samples_sf, 
             file.path(tempdir_r, "samples.shp"), 
             append = FALSE
             )

# Step 4.5 -- Get ground truth points
valid_csv <- system.file(
  "class/samples_validation.csv", package = "sitsdata"
  )

# Step 4.6 -- Calculate accuracy according to Olofsson's method
area_acc <- sits_accuracy(probs_datacube_class, 
                          validation = valid_csv,
                          multicores = 4
                          )

# Step 4.7 -- Print the area estimated accuracy
area_acc

# Step 4.8 -- Show confusion matrix
area_acc$error_matrix