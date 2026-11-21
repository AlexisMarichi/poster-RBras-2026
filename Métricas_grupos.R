############
# Pipeline # 
############

#########################################
# Avaliação quantitativa dos embeddings #
#########################################

# Configuração estrita do Python (Para PHATE e Métricas do Scikit-Learn)
path_to_python <- "/home/users/a291040/miniforge3/envs/sc_rna_seq/bin/python"
Sys.setenv(RETICULATE_PYTHON = path_to_python)
library(reticulate)
use_python(path_to_python, required = TRUE)

# Inicializar e verificar módulos do Python necessários
py_config()
sklearn_manifold <- import("sklearn.manifold")
sklearn_metrics <- import("sklearn.metrics")

####################
# Lista de pacotes #
####################

suppressPackageStartupMessages({
  library(scater) 
  library(SingleCellExperiment) 
  library(scran) 
  library(mclust) # Para ARI
})

#####################
# Reprodutibilidade #
#####################

set.seed(123)

#############################
# Carregar objeto processado #
#############################

sim <- readRDS("splatter_sim_10k_DR.rds")

#####################
# Rótulos verdadeiros #
#####################

true_labels <- colData(sim)$Group
message("Distribuição dos grupos reais (Ground Truth):")
print(table(true_labels))

###################
# Seleção de HVGs #
###################

# Recuperamos os mesmos 2000 HVGs para garantir consistência espacial
dec <- modelGeneVar(sim)
top_hvgs <- getTopHVGs(dec, n = 2000)

####################################################################
# Espaço original de alta dimensão para Trustworthiness            #
####################################################################

message("Preparando espaço original de alta dimensão (2000 HVGs)...")
original_space <- t(as.matrix(logcounts(sim)[top_hvgs, ]))

gc() # Limpeza inicial de memória

####################################
# Função de avaliação otimizada    #
####################################

evaluate_embedding <- function(
    embedding,
    true_labels,
    method_name,
    original_space,
    sklearn_metrics,
    sklearn_manifold
) {
  
  # 1. K-means clustering (Identificação de grupos discretos)
  message(paste("[", method_name, "] Executando K-means..."))
  k <- length(unique(true_labels))
  km <- kmeans(
    embedding,
    centers = k,
    nstart = 20
  )
  pred_labels <- km$cluster
  
  # 2. Adjusted Rand Index (ARI)
  message(paste("[", method_name, "] Calculando Adjusted Rand Index (ARI)..."))
  ari <- mclust::adjustedRandIndex(
    true_labels,
    pred_labels
  )
  
  # 3. Silhouette Score Otimizado (Via Python Scikit-Learn)
  # Evita o cálculo explícito de distâncias pesadas no R
  message(paste("[", method_name, "] Calculando Silhouette Score..."))
  sil_mean <- sklearn_metrics$silhouette_score(
    embedding, 
    as.integer(pred_labels)
  )
  
  # 4. Trustworthiness Otimizado (Via Python Scikit-Learn)
  # Algoritmo indexado em C++ ultrarrápido para vizinhanças locais
  message(paste("[", method_name, "] Calculando Trustworthiness..."))
  trust <- sklearn_manifold$trustworthiness(
    original_space, 
    embedding, 
    n_neighbors = as.integer(15)
  )
  
  gc() # Liberar resíduos temporários de memória do R
  
  # Consolidação da linha
  data.frame(
    Method = method_name,
    ARI = ari,
    Silhouette = sil_mean,
    Trustworthiness = trust
  )
}

########################################
# Execução do Benchmarking             #
########################################

# Avaliação PCA
message("Iniciando avaliação do PCA...")
res_pca <- evaluate_embedding(
  embedding = reducedDim(sim, "PCA"),
  true_labels = true_labels,
  method_name = "PCA",
  original_space = original_space,
  sklearn_metrics = sklearn_metrics,
  sklearn_manifold = sklearn_manifold
)

# Avaliação UMAP
message("Iniciando avaliação do UMAP...")
res_umap <- evaluate_embedding(
  embedding = reducedDim(sim, "UMAP"),
  true_labels = true_labels,
  method_name = "UMAP",
  original_space = original_space,
  sklearn_metrics = sklearn_metrics,
  sklearn_manifold = sklearn_manifold
)

# Avaliação PHATE
message("Iniciando avaliação do PHATE...")
res_phate <- evaluate_embedding(
  embedding = reducedDim(sim, "PHATE"),
  true_labels = true_labels,
  method_name = "PHATE",
  original_space = original_space,
  sklearn_metrics = sklearn_metrics,
  sklearn_manifold = sklearn_manifold
)

########################################
# Combinar e Salvar Resultados         #
########################################

results <- rbind(res_pca, res_umap, res_phate)

message("\n=== TABELA DE RESULTADOS FINAL ===")
print(results)

message("Salvando arquivos de métricas...")
write.csv(results, "metricas_embeddings_10.csv", row.names = FALSE)
saveRDS(results, "metricas_embeddings_10.rds")

########################################
# Encerramento Seguro                  #
########################################

sessionInfo()
gc()
message("--- Avaliação quantitativa finalizada com sucesso total! ---")