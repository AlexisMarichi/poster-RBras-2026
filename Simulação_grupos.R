############
# Pipeline #
############

#######################################
# Estimativa e Simulação com Splatter #
#######################################

###################
#Lista de pacotes #
###################

library(Matrix)
library(SingleCellExperiment)
library(splatter)

##############################################
# Configuração manual do tamanho das células #
##############################################

n_cells_sim <- 10000 # Variable

# Carregar matriz com lógica de cache RDS para economizar tempo

if (file.exists("counts_sub.rds")) {
  message("Carregando matriz real a partir do cache do RDS...")
  counts <- readRDS("counts_sub.rds")
} else {
  message("Lendo arquivos originais MTX...")
  counts <- readMM("counts_sub.mtx")
  genes <- read.delim("genes.tsv", header = FALSE)
  cells <- read.delim("barcodes.tsv", header = FALSE)
  
  # Atribuindo nomes
  rownames(counts) <- genes$V1
  colnames(counts) <- cells$V1
  
  # Converter explicitamente de dgTMatrix para dgCMatrix para que o Splatter não falhe
  counts <- as(counts, "CsparseMatrix")
  
  saveRDS(counts, "counts_sub.rds")
}

########################
# FILTRAR GENES VACÍOS #
########################

# Eliminar genes sin expresión
message("Filtrando genes sin expresión...")
keep_genes <- rowSums(counts) > 0
counts <- counts[keep_genes, ]
gc()

#######################################################################
# Criar SingleCellExperiment (mantendo o formato esparso 'dgCMatrix') #
#######################################################################

# Criamos o objeto passando a lista com o nome exato 'counts'
sce <- SingleCellExperiment(assays = list(counts = counts))

# Atribuímos explicitamente os nomes de genes e células ao objeto SCE.
rownames(sce) <- rownames(counts)
colnames(sce) <- colnames(counts)

# Liberação imediata de RAM
rm(counts) 
gc()

############################
# Estimativa de parâmetros #
############################

message("Estimando parámetros de Splatter...")

# Extraímos a matriz do objeto SCE para passá-la diretamente
# O Splatter aceita uma matriz de genes x células
matriz_para_estimar <- counts(sce) 

# Lançamos a estimativa diretamente sobre a matriz, não sobre o objeto SCE
params <- splatEstimate(as.matrix(matriz_para_estimar)) 
gc()

############################################
# Simulação dos 5 cenários macrobiológicos #  
############################################

set.seed(123) # Para reprodutibilidade

sim <- splatSimulate(
  params,
  batchCells = n_cells_sim,
  method = "groups",
  group.prob = c(0.2, 0.2, 0.2, 0.2, 0.2), # 5 grupos com peso estatístico igual
  verbose = TRUE
)

# Armazenamento dinâmico do resultado
nombre_archivo <- paste0("splatter_sim_", (n_cells_sim / 1000), "k.rds")
saveRDS(sim, nombre_archivo)