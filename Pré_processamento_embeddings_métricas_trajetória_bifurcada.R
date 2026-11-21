#####################################################################
# Avaliação final - Trajetória bifurcada                            #
# PCA vs UMAP vs PHATE                                              #
#####################################################################

########################################
# 1. Pacotes
########################################

suppressPackageStartupMessages({
  
  library(scater)
  library(scran)
  library(phateR)
  library(SingleCellExperiment)
  library(ggplot2)
  library(viridis)
  library(patchwork)
  library(cluster)
  library(mclust)
  library(reticulate)
  
})

set.seed(123)


########################################
# 2. Configuração do Python
########################################

path_to_python <-
  "/home/users/a291040/miniforge3/envs/sc_rna_seq/bin/python"

Sys.setenv(
  RETICULATE_PYTHON = path_to_python
)

use_python(
  path_to_python,
  required = TRUE
)

sklearn_manifold <- import(
  "sklearn.manifold"
)


########################################
# 3. Configuração fixa
########################################

n_hvgs <- 500
n_pcs  <- 6

umap_neighbors <- 15
umap_min_dist  <- 0.01
umap_metric    <- "euclidean"

phate_knn  <- 60
phate_t   <- 22
phate_seed <- 123


########################################
# 4. Carregar simulação
########################################

archivo_base <-
  "splatter_sim_10k_bifurc_path.rds"

if (!file.exists(archivo_base)) {
  
  stop(
    paste(
      "Arquivo não encontrado:",
      archivo_base
    )
  )
  
}

message(
  "\nCarregando simulação bifurcada..."
)

sim <- readRDS(
  archivo_base
)


########################################
# 5. Inspeção da estrutura
########################################

message(
  "\nDistribuição dos paths:"
)

print(
  table(colData(sim)$Group)
)


########################################
# 6. Normalização
########################################

message(
  "\nCalculando log-normalização..."
)

sim <- logNormCounts(sim)


########################################
# 7. Seleção de HVGs
########################################

message(
  paste(
    "\nSelecionando os",
    n_hvgs,
    "HVGs..."
  )
)

dec <- modelGeneVar(sim)

top_hvgs <- getTopHVGs(
  dec,
  n = n_hvgs
)


########################################
# 8. Espaço original para
#    Trustworthiness
########################################

original_space <-
  t(
    as.matrix(
      logcounts(sim)[top_hvgs, ]
    )
  )


########################################
# 9. PCA
########################################

message(
  "\nExecutando PCA..."
)

sim <- runPCA(
  sim,
  subset_row = top_hvgs,
  ncomponents = n_pcs
)


########################################
# 10. UMAP
########################################

message(
  "\nExecutando UMAP..."
)

sim <- runUMAP(
  sim,
  dimred = "PCA",
  n_neighbors = umap_neighbors,
  min_dist = umap_min_dist,
  metric = umap_metric,
  seed = 123
)


########################################
# 11. PHATE
########################################

message(
  "\nExecutando PHATE..."
)

n_cores <-
  as.integer(
    Sys.getenv(
      "SLURM_CPUS_PER_TASK",
      unset = 1
    )
  )

phate_out <- phate(
  reducedDim(sim, "PCA"),
  knn = phate_knn,
  t = phate_t,
  n.jobs = n_cores,
  seed = phate_seed
)

reducedDim(sim, "PHATE") <-
  as.matrix(phate_out)

rm(phate_out)

gc()


######################################################################
# 12. Definir os ramos verdadeiros
######################################################################

message(
  "\nDefinindo os ramos verdadeiros..."
)

groups <- as.character(
  colData(sim)$Group
)

idx_branch <- which(
  groups %in% c(
    "Path2",
    "Path3"
  )
)

branch_true <- factor(
  groups[idx_branch],
  levels = c(
    "Path2",
    "Path3"
  )
)

message(
  "\nNúmero de células por ramo:"
)

print(
  table(branch_true)
)


######################################################################
# 13. Construção da distância verdadeira
#     sobre a árvore bifurcada
######################################################################

message(
  "\nConstruindo distâncias verdadeiras sobre a árvore..."
)

steps <- as.numeric(
  colData(sim)$Step
)

groups_all <- as.character(
  colData(sim)$Group
)


# Nó de bifurcação:
#
# Path1 termina em Step = 100
# Path2 e Path3 começam a partir desse nó.

branch_point <- 100


# Distância geodésica entre duas células
# dentro da árvore.

tree_distance <- function(
    group_i,
    step_i,
    group_j,
    step_j
) {
  
  # Mesmo ramo
  if (group_i == group_j) {
    
    return(
      abs(step_i - step_j)
    )
    
  }
  
  
  # Path1 -> Path2
  if (
    group_i == "Path1" &&
    group_j == "Path2"
  ) {
    
    return(
      (branch_point - step_i) +
        step_j
    )
    
  }
  
  
  # Path2 -> Path1
  if (
    group_i == "Path2" &&
    group_j == "Path1"
  ) {
    
    return(
      step_i +
        (branch_point - step_j)
    )
    
  }
  
  
  # Path1 -> Path3
  if (
    group_i == "Path1" &&
    group_j == "Path3"
  ) {
    
    return(
      (branch_point - step_i) +
        step_j
    )
    
  }
  
  
  # Path3 -> Path1
  if (
    group_i == "Path3" &&
    group_j == "Path1"
  ) {
    
    return(
      step_i +
        (branch_point - step_j)
    )
    
  }
  
  
  # Path2 -> Path3
  if (
    group_i == "Path2" &&
    group_j == "Path3"
  ) {
    
    return(
      step_i +
        step_j
    )
    
  }
  
  
  # Path3 -> Path2
  if (
    group_i == "Path3" &&
    group_j == "Path2"
  ) {
    
    return(
      step_i +
        step_j
    )
    
  }
  
  return(NA)
  
}


########################################
# 14. Amostra para distância global
########################################

set.seed(123)

n_cells <- ncol(sim)

sample_size <-
  min(
    1000,
    n_cells
  )

idx_sample <-
  sample(
    seq_len(n_cells),
    size = sample_size
  )


########################################
# Matriz verdadeira
########################################

dist_true <- matrix(
  0,
  nrow = sample_size,
  ncol = sample_size
)

for (i in seq_len(sample_size)) {
  
  for (j in seq_len(sample_size)) {
    
    dist_true[i, j] <-
      tree_distance(
        groups_all[idx_sample[i]],
        steps[idx_sample[i]],
        groups_all[idx_sample[j]],
        steps[idx_sample[j]]
      )
    
  }
  
}


######################################################################
# 15. Função de avaliação
######################################################################

evaluate_embedding <- function(
    embedding,
    method_name
) {
  
  
  ################################################
  # 1. Branch ARI
  ################################################
  
  embedding_branch <-
    embedding[idx_branch, , drop = FALSE]
  
  
  set.seed(123)
  
  km <- kmeans(
    embedding_branch,
    centers = 2,
    nstart = 50
  )
  
  
  branch_ari <-
    adjustedRandIndex(
      km$cluster,
      branch_true
    )
  
  
  ################################################
  # 2. Distância global
  ################################################
  
  embedding_sample <-
    embedding[idx_sample, , drop = FALSE]
  
  
  dist_embedding <-
    as.matrix(
      dist(
        embedding_sample
      )
    )
  
  
  rho_global_dist <-
    cor(
      as.vector(dist_true),
      as.vector(dist_embedding),
      method = "spearman"
    )
  
  
  ################################################
  # 3. Trustworthiness
  ################################################
  
  trust <-
    sklearn_manifold$trustworthiness(
      original_space,
      embedding,
      n_neighbors = as.integer(15)
    )
  
  
  ################################################
  # Resultado
  ################################################
  
  data.frame(
    
    HVG_Count = n_hvgs,
    
    PC_Count = n_pcs,
    
    Method = method_name,
    
    Branch_ARI =
      round(
        branch_ari,
        3
      ),
    
    Global_Distance_Spearman =
      round(
        rho_global_dist,
        3
      ),
    
    Trustworthiness =
      round(
        trust,
        3
      )
    
  )
  
}


######################################################################
# 16. Avaliação dos três métodos
######################################################################

message(
  "\nCalculando métricas..."
)


res_pca <-
  evaluate_embedding(
    reducedDim(
      sim,
      "PCA"
    )[, 1:2],
    "PCA"
  )


res_umap <-
  evaluate_embedding(
    reducedDim(
      sim,
      "UMAP"
    ),
    "UMAP"
  )


res_phate <-
  evaluate_embedding(
    reducedDim(
      sim,
      "PHATE"
    ),
    "PHATE"
  )


results <-
  rbind(
    res_pca,
    res_umap,
    res_phate
  )


######################################################################
# 17. Exibir resultados
######################################################################

message(
  "\n=================================================="
)

message(
  "RESULTADOS - TRAJETÓRIA BIFURCADA"
)

message(
  "=================================================="
)

print(
  results
)


######################################################################
# 18. Salvar métricas
######################################################################

nombre_csv <-
  "metricas_bifurcacion_HVG_500_PC_6.csv"

write.csv(
  results,
  nombre_csv,
  row.names = FALSE
)


######################################################################
# 19. Visualização
######################################################################

message(
  "\nGerando visualizações..."
)


var_color <- "Group"


p1 <-
  plotReducedDim(
    sim,
    "PCA",
    colour_by = var_color,
    point_size = 1
  ) +
  ggtitle(
    "PCA"
  ) +
  theme_classic()


p2 <-
  plotReducedDim(
    sim,
    "UMAP",
    colour_by = var_color,
    point_size = 1
  ) +
  ggtitle(
    "UMAP"
  ) +
  theme_classic()


p3 <-
  plotReducedDim(
    sim,
    "PHATE",
    colour_by = var_color,
    point_size = 1
  ) +
  ggtitle(
    "PHATE"
  ) +
  theme_classic()


######################################################################
# 20. Salvar PDF
######################################################################

nombre_pdf <-
  "Visualizacao_final_trajetoria_bifurcada.pdf"


pdf(
  nombre_pdf,
  width = 16,
  height = 5
)


print(
  p1 + p2 + p3 +
    plot_layout(
      guides = "collect"
    ) &
    theme(
      legend.position = "right"
    )
)


invisible(
  dev.off()
)


######################################################################
# 21. Salvar objeto
######################################################################

nombre_rds <-
  "splatter_sim_bifurc_DR_HVG_500_PC_6.rds"


saveRDS(
  sim,
  nombre_rds
)


######################################################################
# FIM
######################################################################

message(
  "\n=================================================="
)

message(
  "ANÁLISE DA TRAJETÓRIA BIFURCADA FINALIZADA"
)

message(
  "=================================================="
)

message(
  paste(
    "Métricas:",
    nombre_csv
  )
)

message(
  paste(
    "Figura:",
    nombre_pdf
  )
)

message(
  paste(
    "Objeto:",
    nombre_rds
  )
)

message(
  "=================================================="
)