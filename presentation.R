# ------------------------------------------------------------
# 1. Instalação e carregamento dos pacotes
# ------------------------------------------------------------

library(sits)       # Pacote principal do Brazil Data Cube para análise e classificação de séries temporais de imagens de satélite.
library(sitsdata)   # Pacote auxiliar que fornece exemplos prontos de cubos de dados e séries temporais (útil para testes e aprendizado).
library(sf)         # Manipula dados vetoriais espaciais (como shapefiles e GeoJSON), permitindo operações espaciais modernas.
library(tibble)     # Estrutura moderna de tabelas, mais legível e segura que data.frame.
library(dplyr)      # Facilita operações com dados tabulares (como filtrar, agrupar, selecionar e modificar colunas).
library(rstac)      # Permite acessar catálogos de dados no padrão STAC, como o Brazil Data Cube (BDC).

# ------------------------------------------------------------
# 2. Criação do cubo de dados (Sentinel-2)
# ------------------------------------------------------------

# Cria o cubo de dados Sentinel-2 a partir do provedor "BDC" (Brazil Data Cube).
# O parâmetro "collection" define o produto: SENTINEL-2-16D é o mosaico de 16 dias.
# As bandas selecionadas incluem bandas espectrais e índices de vegetação (NDVI, EVI) e máscara de nuvens (CLOUD).
# O intervalo temporal vai de 2020 a 2024.
cube <- sits_cube(
  source = "BDC",                               # Define o provedor como o Brazil Data Cube.
  collection = "SENTINEL-2-16D",                # Especifica o produto Sentinel-2 com composição de 16 dias.
  bands = c("B02", "B03", "B04", "B08", "B8A", "B11", "NDVI", "EVI", "CLOUD"), # Seleciona bandas relevantes para classificação.
  tiles = c("012014"),                          # Define o tile da área de interesse (Rondônia neste caso).
  start_date = "2020-01-01",                    # Data inicial das imagens.
  end_date = "2024-12-31",                      # Data final das imagens.
  progress = TRUE                               # Exibe barra de progresso durante o download e processamento.
)

# Filtra novamente o cubo, garantindo que usamos apenas as bandas e o intervalo desejado.
cube_select <- sits_select(
  cube,
  bands = c("B02", "B03", "B04", "B08", "B8A", "B11", "NDVI", "EVI", "CLOUD"), # Bandas selecionadas para análise.
  tiles = c("012014"),                          # Mesmo tile.
  start_date = "2020-01-01",                    # Mesmo intervalo temporal.
  end_date = "2024-12-31",
  progress = TRUE
)

# Exibe uma composição colorida para inspeção visual.
# B11 (SWIR), B8A (NIR) e B02 (Blue) ajudam a distinguir vegetação, solo exposto e corpos d'água.
plot(cube_select, red = "B11", green = "B8A", blue = "B02")

# ------------------------------------------------------------
# 3. Leitura e análise das amostras
# ------------------------------------------------------------

# Lê o shapefile que contém as amostras de treinamento (cada ponto/polígono possui uma classe associada).
samples_sf <- st_read(samp_path)

# Exibe o número total de amostras no arquivo shapefile.
cat("Total de amostras:", nrow(samples_sf), "\n")

# Mostra as dimensões do arquivo (linhas x colunas).
cat("Dimensões (linhas x colunas):", dim(samples_sf), "\n")

# Exibe a quantidade de amostras por classe original.
table(samples_sf$label)

# Padroniza os nomes das classes para uma nomenclatura mais clara e internacional.
samples_sf <- samples_sf %>%
  mutate(label = recode(
    label,
    "dms_queimada"      = "Clear Cut Burned Area",  # Área recém-desmatada e queimada.
    "dms_veg"           = "Clear Cut Vegetation",   # Vegetação em regeneração inicial após desmate.
    "dms_solo_exposto"  = "Clear Cut Bare Soil",    # Solo exposto após desmate.
    "floresta"          = "Forest",                 # Floresta primária.
    "hidrografia"       = "Water bodies",           # Corpos d’água.
    "floresta_alagavel" = "Wetland",                # Floresta alagável.
    "nao_floresta"      = "Non-forest"              # Áreas não florestais (agricultura, pastagem, etc.).
  ))

# Visualiza as amostras no mapa para verificar sua distribuição espacial.
sits_view(samples_sf)

# ------------------------------------------------------------
# 4. Extração das Séries Temporais
# ------------------------------------------------------------

# Extrai séries temporais de reflectância e índices (NDVI, EVI) para cada amostra a partir do cubo Sentinel-2.
# "multicores = 12" usa 12 núcleos do processador para acelerar o processamento.
samples_rondonia_test <- sits_get_data(
  cube = cube_select,                # Cubo de dados de entrada.
  samples = samples_sf,              # Conjunto de amostras com classes conhecidas.
  start_date = "2020-01-01",         # Data inicial da série temporal.
  end_date = "2024-12-31",           # Data final da série temporal.
  label = "label",                   # Coluna com os rótulos das classes.
  multicores = 12,                   # Número de núcleos de CPU usados no processamento paralelo.
  progress = TRUE                    # Mostra progresso.
)

# Seleciona apenas as bandas NDVI e EVI e plota os padrões médios temporais de cada classe.
samples_rondonia_test |>
  sits_select(bands = c("NDVI", "EVI"), start_date = "2020-01-01", end_date = "2024-12-18") |>  # Seleciona bandas e período.
  sits_patterns() |>  # Gera os padrões médios das séries temporais por classe.
  plot()              # Plota os padrões (útil para verificar separabilidade entre classes).

# ------------------------------------------------------------
# 5. Avaliação da Qualidade das Amostras (SOM)
# ------------------------------------------------------------

# Cria um mapa SOM (Self-Organizing Map) que organiza as séries temporais em grupos similares.
# Essa etapa ajuda a identificar amostras inconsistentes ou com sobreposição entre classes.
som_cluster <- sits_som_map(
  samples_rondonia_test,   # Séries temporais das amostras.
  grid_xdim = 12,          # Número de células no eixo X do grid SOM (define resolução da análise).
  grid_ydim = 12,          # Número de células no eixo Y do grid SOM.
  rlen = 100,              # Número de iterações de treinamento do SOM.
  distance = "dtw",        # Usa distância DTW (Dynamic Time Warping) para medir similaridade temporal.
  som_radius = 2,          # Define o raio de vizinhança (influência entre células).
  mode = "online"          # Modo de atualização do SOM (mais rápido).
)

# Exibe o mapa SOM, onde cores semelhantes indicam séries temporais parecidas.
plot(som_cluster)

# Avalia a pureza dos clusters, identificando classes misturadas ou inconsistentes.
som_eval <- sits_som_evaluate_cluster(som_cluster)
som_eval   # Mostra tabela com resultados de avaliação dos clusters.

# ------------------------------------------------------------
# 6. Treinamento e Classificação
# ------------------------------------------------------------

# Define uma semente aleatória para garantir reprodutibilidade dos resultados.
set.seed(03022024)

# Treina um modelo de Random Forest usando as amostras temporais.
# O algoritmo cria múltiplas árvores de decisão e combina seus resultados para reduzir erros.
rf_model <- sits_train(
  samples = samples_rondonia_test,   # Séries temporais com rótulos.
  ml_method = sits_rfor()            # Define o método de aprendizado como Random Forest.
)

# Exibe a importância das variáveis (bandas e índices) no modelo.
plot(rf_model)

# ------------------------------------------------------------
# 7. Classificação do Cubo de Dados
# ------------------------------------------------------------

# Aplica o modelo Random Forest ao cubo para classificar cada pixel.
# O resultado é um mapa de probabilidades por classe.
class_prob <- sits_classify(
  data = cube_select,                # Cubo de entrada (imagens).
  ml_model = rf_model,               # Modelo treinado.
  output_dir = classification_path,  # Pasta onde os resultados serão salvos.
  version = "train",                 # Identificador da versão da classificação.
  multicores = 26,                   # Usa 26 núcleos para processamento paralelo (ajuste conforme hardware).
  memsize = 50,                      # Tamanho máximo de memória (GB) permitido.
  progress = TRUE                    # Exibe progresso na classificação.
)

# Visualiza o mapa de probabilidades por classe.
sits_view(class_prob)

# ------------------------------------------------------------
# 8. Geração do Mapa Final
# ------------------------------------------------------------

# Converte o mapa de probabilidades em um mapa de classes, 
# atribuindo a cada pixel a classe com maior probabilidade.
class_map <- sits_label_classification(
  cube = class_prob,                 # Mapa de probabilidades.
  output_dir = classification_path,  # Diretório de saída.
  version = "train",                 # Versão de classificação.
  multicores = 26,                   # Núcleos usados.
  memsize = 50,                      # Memória alocada (GB).
  progress = TRUE                    # Exibe progresso.
)

# Visualiza o mapa final de classificação com as classes temáticas.
sits_view(class_map)
