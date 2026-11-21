# Pacotes
import scanpy as sc 
from scipy.io import mmwrite
import scipy.sparse as sp
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split

# Conectar ao h5ad no modo "Backed"
adata = sc.read_h5ad("All_Neurons.h5ad", backed="r")

# Definir o tamanho da subamostra para a estimativa
target_size = 30000  
strata_col = "supercluster_term"

# Criar um vetor de índices sequenciais (0 a N-1)
indices = np.arange(adata.n_obs)

# Realizar a partição estratificada
# O argumento 'stratify' garante que as proporções dos superclusters sejam mantidas.
_, idx_sub = train_test_split(
    indices, 
    test_size=target_size, 
    stratify=adata.obs[strata_col], 
    random_state=123
)

# Organizar os índices para que o HDF5 realize a leitura de forma sequencial e ultrarrápida
idx_sub = np.sort(idx_sub)

# Extrair a subamostra virtual
adata_sub = adata[idx_sub, :]
print("\nSubamostra criada com sucesso:")
print(adata_sub)

# Verificar se as proporções foram mantidas
print("\nDistribuição percentual na subamostra (Top 5):")
print(adata_sub.obs[strata_col].value_counts(normalize=True).head())

# Extrair contagens (Priorizando dados brutos, se existirem)
if adata_sub.raw is not None:
    X = adata_sub.raw.X
    genes = adata_sub.raw.var_names
    print("Extraindo dados do slot: adata.raw.X")
else:
    X = adata_sub.X
    genes = adata_sub.var_names
    print("Extraindo dados do slot: adata.X")

cells = adata_sub.obs_names

# Converter a matriz virtual 'backed' em uma matriz esparsa real na RAM
if hasattr(X, "to_memory"):
    X = X.to_memory()
elif not sp.issparse(X):
    X = sp.csr_matrix(X)

# Transpor para o R (Passar de Células x Genes para Genes x Células)
X_to_save = X.T
print(f"Dimensões prontas para o R (Genes): {X_to_save.shape[0]}, Células: {X_to_save.shape[1]})")

# Salvar arquivos
mmwrite("counts_sub.mtx", X_to_save)

pd.Series(genes).to_csv("genes.tsv", sep="\t", index=False, header=False)
pd.Series(cells).to_csv("barcodes.tsv", sep="\t", index=False, header=False)

# Salve também os metadados das células selecionadas, caso precise deles no R
adata_sub.obs.to_csv("metadata_sub.tsv", sep="\t")
