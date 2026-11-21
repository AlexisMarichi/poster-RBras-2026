#####################################################################
# Pipeline Unificado: Redução de Dimensionalidade + Métricas HVG    #
#####################################################################

# Capturar número de HVGs a partir da linha de comando
args <- commandArgs(trailingOnly = TRUE)
n_hvgs <- if (length(args) > 0) as.integer(args[1]) else 1000

message(paste0("\n=================================================="))
message(paste0("   INICIANDO PIPELINE COM: ", n_hvgs, " HVGs"))
message(paste0("=================================================="))

######################################
# Configuração e Importação Python   #
######################################

path_to_python <- "/home/users/a291040/miniforge3/envs/sc_rna_seq/bin/python"
Sys.setenv(RETICULATE_PYTHON = path_to_python)
library(reticulate)
use_python(path_to_python, required = TRUE)

sklearn_manifold <- import("sklearn.manifold")

####################
# Lista de Pacotes #
####################

suppressPackageStartupMessages({
  library(scater) 
  library(scran) 
  library(phateR) 
  library(SingleCellExperiment) 
  library(patchwork)  
  library(ggplot2)
  library(viridis)
  library(princurve)
})

set.seed(123)

############################
# Carregar Objeto Simulado #
############################

archivo_base <- "splatter_sim_10k_linear_path.rds"
if (!file.exists(archivo_base)) {
  stop(paste("Erro! Arquivo não encontrado:", archivo_base))
}

sim <- readRDS(archivo_base)

# Reconstrução matemática do Pseudotempo Contínuo (1 a 400)
raw_path <- as.numeric(gsub("[^0-9]", "", as.character(colData(sim)$Group)))
raw_step <- as.numeric(gsub("[^0-9]", "", as.character(colData(sim)$Step)))

# Cálculo do tempo linear contínuo acumulado
true_step <- (raw_path - 1) * 100 + raw_step

# Salvar o tempo global para avaliação e visualização
colData(sim)$Global_Step <- true_step

# true_step serve para realizar os cálculos matemáticos diretamente no código 
# (como o teste de Spearman)

####################################
# Normalização e Seleção de HVGs   #
####################################

message("Calculando Log-Normalização...")
sim <- logNormCounts(sim)

message(paste("Selecionando top", n_hvgs, "HVGs..."))
dec <- modelGeneVar(sim)
top_hvgs <- getTopHVGs(dec, n = n_hvgs)

# Matriz para o cálculo de Trustworthiness (Espaço Original)
original_space <- t(as.matrix(logcounts(sim)[top_hvgs, ]))

# Subamostra para Preservação de Distâncias (Eficiência de Memória) 
set.seed(123)
n_cells <- ncol(sim)
# Garante selecionar no máximo 1000 células aleatórias para não estourar a RAM em datasets grandes
sample_size <- min(1000, n_cells)
idx_sample <- sample(seq_len(n_cells), size = sample_size)

true_step_sub <- true_step[idx_sample]
dist_true_step <- as.matrix(dist(true_step_sub)) # Matriz de distâncias em Step

########################################
# Redução de Dimensionalidade (DR)     #
########################################

# PCA
message("Executando PCA (10 PCs)...")
sim <- runPCA(sim, subset_row = top_hvgs, ncomponents = 10)

# UMAP
message("Executando UMAP...")
sim <- runUMAP(sim, dimred = "PCA", n_neighbors = 15, min_dist = 0.01, metric = "euclidean")

# PHATE
message("Executando PHATE...")
phate_out <- phate(reducedDim(sim, "PCA"), n.jobs = n_cores, seed = 123, knn = 5)
reducedDim(sim, "PHATE") <- as.matrix(phate_out)
rm(phate_out)
gc()

########################################
# Geração de Visualização Conjunta     #
########################################

message("Gerando gráficos comparativos...")
# var_color <- "Step"
var_color <- "Global_Step" 

p1 <- plotReducedDim(sim, "PCA", colour_by = var_color, point_size = 0.5) + 
  scale_color_viridis_c(option = "magma") + ggtitle(paste("PCA (", n_hvgs, "HVGs)")) + theme_classic()

p2 <- plotReducedDim(sim, "UMAP", colour_by = var_color, point_size = 0.5) + 
  scale_color_viridis_c(option = "magma") + ggtitle(paste("UMAP (", n_hvgs, "HVGs)")) + theme_classic()

p3 <- plotReducedDim(sim, "PHATE", colour_by = var_color, point_size = 0.5) + 
  scale_color_viridis_c(option = "magma") + ggtitle(paste("PHATE (", n_hvgs, "HVGs)")) + theme_classic()

pdf(paste0("Visualizacion_conjunta_HVG_", n_hvgs, ".pdf"), width = 16, height = 5)
print((p1 + p2 + p3) + plot_layout(guides = "collect") & theme(legend.position = "right"))
dev.off()

########################################
# Avaliação Quantitativa de Métricas   #
########################################

evaluate_embedding <- function(embedding, method_name) {
  
  # Correlação de Pseudotempo 1D (Curva Principal)
  
  # A) Criamos uma matriz guia inicial ('start') tomando ~50 pontos ordenados pelo tempo real
  # Isso ensina ao princurve a geometria real (o arco/parábola) para que não corte caminho pelo centro
  ord_idx <- order(true_step)
  step_seq <- round(seq(1, length(ord_idx), length.out = 50))
  start_matrix <- as.matrix(embedding[ord_idx[step_seq], ])
  
  # B) Ajustamos a curva principal usando a matriz guia
  fit_curve <- princurve::principal_curve(
    as.matrix(embedding), 
    start = start_matrix
  )
  
  inferred_pseudotime <- fit_curve$lambda
  
  # C) Orientação da raiz: Garantimos que o pseudotempo comece próximo de 0 na raiz real
  root_idx <- which.min(true_step)
  if (inferred_pseudotime[root_idx] > max(inferred_pseudotime) / 2) {
    inferred_pseudotime <- max(inferred_pseudotime) - inferred_pseudotime
  }
  
  # D) Calculamos a correlação de Spearman
  rho_pseudotime <- abs(cor(inferred_pseudotime, true_step, method = "spearman"))
  
  # Preservação de Distâncias Globais
  dist_embedding <- as.matrix(dist(embedding[idx_sample, ]))
  rho_global_dist <- cor(as.vector(dist_true_step), as.vector(dist_embedding), method = "spearman")
  
  # Trustworthiness
  trust <- sklearn_manifold$trustworthiness(original_space, embedding, n_neighbors = as.integer(15))
  
  data.frame(
    HVG_Count = n_hvgs,
    Method = method_name,
    Pseudotime_Spearman = round(rho_pseudotime, 3),
    Global_Distance_Spearman = round(rho_global_dist, 3),
    Trustworthiness = round(trust, 3)
  )
}

message("Calculando métricas quantitativas...")

#Avaliação justa do PCA em 2D
res_pca   <- evaluate_embedding(reducedDim(sim, "PCA")[, 1:2], "PCA")
res_umap  <- evaluate_embedding(reducedDim(sim, "UMAP"), "UMAP")
res_phate <- evaluate_embedding(reducedDim(sim, "PHATE"), "PHATE")

results <- rbind(res_pca, res_umap, res_phate)

# Salvar Resultados Quantitativos
write.csv(results, paste0("metricas_HVG_", n_hvgs, ".csv"), row.names = FALSE)

# Salvar Objeto Processado caso seja necessária análise posterior
saveRDS(sim, paste0("splatter_sim_linear_DR_HVG_", n_hvgs, ".rds"))

message(paste0("=================================================="))
message(paste0("   SUCESSO: Processamento finalizado para HVG: ", n_hvgs))
message(paste0("=================================================="))