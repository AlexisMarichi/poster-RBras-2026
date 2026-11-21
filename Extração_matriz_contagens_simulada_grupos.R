library(SingleCellExperiment)
library(Matrix)

# Carregar o objeto simulado que o seu pipeline gerou
sim <- readRDS("splatter_sim_10k.rds")

# Extrair a matriz de contagens brutas (mantém o formato eficiente dgCMatrix)
matriz_conteos_sim <- counts(sim)

# Salvar diretamente como uma matriz independente
saveRDS(matriz_conteos_sim, "matriz_conteos_original_simulada_10k.rds")

# Controle de qualidade para garantir que tudo esteja certo antes do envio
print(dim(matriz_conteos_sim)) # Deve mostrar: 42206 genes 10000 células