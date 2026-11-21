############
# Pipeline #
############

##################################
# Normalização e seleção de HVGs #
##################################

# Configuração estrita do Python (será ativada apenas ao executar o PHATE)
path_to_python <- "/home/users/a291040/miniforge3/envs/sc_rna_seq/bin/python"
Sys.setenv(RETICULATE_PYTHON = path_to_python)
library(reticulate)
use_python(path_to_python, required = TRUE)
py_config()

####################
# Lista de pacotes #
####################

suppressPackageStartupMessages({
  library(scater) 
  library(scran) 
  library(phateR) 
  library(SingleCellExperiment) 
  library(patchwork)  
  library(ggplot2)
  library(Matrix)
})

#####################
# Reprodutibilidade #
#####################

set.seed(123)

#######################
# Detecção de núcleos #
#######################

# Detectar núcleos (cores) alocados pelo SLURM para acelerar o PHATE
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))
message(paste("Usando", n_cores, "núcleos para os cálculos paralelos..."))

############################
# Carregar objeto simulado #
############################

sim <- readRDS("splatter_sim_10k.rds")

####################
# Verificar objeto #
####################

print(dim(sim))
print(assayNames(sim))
print(table(colData(sim)$Group)) # Os rótulos verdadeiros estão em colData(sim)$Group

gc() # Libera a memória RAM ocupada por objetos que não estão mais em uso

#####################
# Log-Normalização  #
#####################

message("Calculando Log-Normalização...")
sim <- logNormCounts(sim) # Normalizado com scater

gc()

###################
# Seleção de HVGs #
###################

message("Selecionando HVGs...")

#Primeiro critério 
dec <- modelGeneVar(sim) # Seleção de HVGs via scran
top_hvgs <- getTopHVGs(dec, n = 2000)
message(paste("HVGs selecionados:", length(top_hvgs)))

#Segundo critério (Descomentar este código para usar o critério de seleção FDR de HVGs)
# dec <- modelGeneVar(sim) 
# hvg_indices <- which(dec$FDR < 0.05)
# top_hvgs <- rownames(dec)[hvg_indices]
# message(paste("HVGs seleccionados:", length(top_hvgs)))

#######
# PCA #
#######

message("Executando PCA...")

sim <- runPCA(sim, subset_row = top_hvgs, ncomponents = 30)

#Primer critério
# Gráfico de Codo (Elbow Plot)
message("Generando Elbow Plot...")
# Extraer el porcentaje de varianza explicada por cada PC
percent_var <- attr(reducedDim(sim, "PCA"), "percentVar")

# Crear un data frame para ggplot
elbow_data <- data.frame(PC = 1:length(percent_var), Varianza = percent_var)

# Graficar el codo
p_elbow <- ggplot(elbow_data, aes(x = PC, y = Varianza)) +
  geom_point(color = "blue", size = 1.5) +
  geom_line(color = "darkgray") +
  scale_y_log10() + 
  labs(title = "Elbow Plot (Escala logarítmica)",
       x = "Componentes Principais (PCs)",
       y = "% de Variância Explicada (Escala Log)") +
  theme_classic()

message("Gerando gráfico do elbow...")
pdf("PCA_elbow_plot_log.pdf", width = 7, height = 5)
print(p_elbow)
dev.off()

gc()

#Segundo critério
#sim <- denoisePCA(sim, technical = dec, subset.row = top_hvgs)

########
# UMAP #
########

message("Executando UMAP...")
sim <- runUMAP(sim, dimred = "PCA")

gc()

#########
# PHATE #
#########

message("Executando PHATE (isso pode levar algum tempo)...")
# Adiciona-se n.jobs para aproveitar o paralelismo do cluster
phate_out <- phate(reducedDim(sim, "PCA"), n.jobs = n_cores)

# Converter embedding para matriz explicitamente
# phate_embedding <- as.matrix(phate_out$embedding)
# OTIMIZAÇÃO: Extração segura da matriz de embedding
phate_embedding <- as.matrix(phate_out)

# Salvar embedding PHATE no objeto para manter consistência com runPCA/runUMAP
reducedDim(sim, "PHATE") <- phate_embedding

# Liberar memória
rm(phate_out, phate_embedding)
gc()

####################
# Visualização PCA #
####################

pdf("PCA_groups.pdf", width = 8, height = 6)
print(plotReducedDim(sim, dimred = "PCA", colour_by = "Group", point_size = 0.3))
dev.off()

####################
# Visualização UMAP #
####################

pdf("UMAP_groups.pdf", width = 8, height = 6)
print(plotReducedDim(sim, dimred = "UMAP", colour_by = "Group", point_size = 0.3))
dev.off()

#####################
# Visualização PHATE #
#####################

pdf("PHATE_groups.pdf", width = 8, height = 6)
print(plotReducedDim(sim, dimred = "PHATE", colour_by = "Group", point_size = 0.3))
dev.off()

#########################
# Visualização conjunta #
#########################

message("Gerando gráficos conjuntos...")

p1 <- plotReducedDim(sim, "PCA", colour_by = "Group", point_size = 0.3) + ggtitle("PCA") + theme_classic()
p2 <- plotReducedDim(sim, "UMAP", colour_by = "Group", point_size = 0.3) + ggtitle("UMAP") + theme_classic()
p3 <- plotReducedDim(sim, "PHATE", colour_by = "Group", point_size = 0.3) + ggtitle("PHATE") + theme_classic()

pdf("Visualizacion_conjunta.pdf", width = 15, height = 5)
print((p1 + p2 + p3) + plot_layout(guides = "collect") & theme(legend.position = "right"))
dev.off()

#######################
# Salvar objeto final #
#######################

message("Salvando objeto final...")
saveRDS(sim, "splatter_sim_10k_DR.rds") # Conterá Normalized Counts, PCA, UMAP e PHATE

######################
# Informações finais #
######################

sessionInfo()
gc()
message("--- Pipeline finalizado com sucesso! ---")