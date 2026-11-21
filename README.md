<div align="justify">
    
# poster-RBras-2026
Dados e códigos do pôster comparativo de redução de dimensionalidade em scRNA-seq - 70ª RBras / IMECC / UNICAMP

# Conjunto de Dados de Simulação e Benchmark de Métodos de Redução de Dimensionalidade em scRNA-seq (PCA, UMAP e PHATE)

## Autores e Afiliação
- **Autor:** Alexis Marichi Zela (IMECC - UNICAMP)
- **Orientador:** Prof. Dr. Benilton de Sá Carvalho (IMECC - UNICAMP)
- **Instituição:** Instituto de Matemática, Estatística e Computação Científica (IMECC), Universidade Estadual de Campinas (UNICAMP).

---

## 1. Descrição Geral
Este conjunto de dados e código-fonte acompanham a dissertação de mestrado intitulada **"Métodos Estatísticos de Redução de Dimensionalidade em scRNA-seq: Fundamentos, Avaliação e Recomendações Práticas"** (UNICAMP, 2026). 

O repositório contém as matrizes de contagens simuladas, tabelas de métricas quantitativas e scripts de reprodução dos experimentos de *benchmark* comparativo entre os métodos de redução de dimensionalidade **PCA**, **UMAP** e **PHATE**. As simulações foram geradas utilizando o pacote `Splatter` (modelo `Splat`), calibrado a partir do atlas neuronal humano **"All Neurons"** (Siletti et al., 2023 / Human Cell Atlas).

---

</div>

```text
├── README.md
├── Scripts/
│   ├── Amostragem_estratificada.py                                   # Script Python para amostragem estratificada do atlas .h5ad
│   ├── Simulação_grupos.R                                            # Simulação de populações celulares discretas
│   ├── Simulação_trajetória_linear.R                                 # Simulação de trajetória contínua linear
│   ├── Simulação_trajetória_bifurcada.R                              # Simulação de trajetória contínua bifurcada
│   ├── Visualizar_parâmetros_simulação_grupos.R                      # Análise e visualização dos parâmetros simulados (adaptável para o caso de trajetórias)
│   ├── Extração_matriz_contagens_simulada_grupos.R                   # Extração e preparação das matrizes base (adaptável para o caso de trajetórias)
│   ├── Pré_processamento_embeddings_grupos.R                         # Execução de DR (PCA, UMAP, PHATE) para grupos discretos (adaptável para os outros parâmetros)
│   ├── Pré_processamento_embeddings_métricas_trajetória_linear.R     # Execução de DR e cálculo de métricas para trajetória linear (adaptável para os outros parâmetros)
│   ├── Pré_processamento_embeddings_métricas_trajetória_bifurcada.R  # Execução de DR e cálculo de métricas para bifurcação (validação da configuração ótima)
│   ├── Métricas_grupos.R                                             # Cálculo de métricas (ARI, Silhouette, etc.) para grupos
│   └── Gráficos_tendência_trajetória_linear.R                        # Visualização do impacto de parâmetros (ex: número de HVGs) (adaptável para os outros parâmetros) 
└── Dados_e_Resultados/
    ├── Matriz_contagens_amostrada.rds                                # Dados reais amostrados do atlas de referência
    ├── Matriz_contagens_simulada_grupos.rds                          # Objeto sintético do cenário de grupos discretos
    ├── Matriz_contagens_simulada_trajetória_linear.rds               # Objeto sintético do cenário de trajetória linear
    ├── Matriz_contagens_simulada_trajetória_bifurcada.rds            # Objeto sintético do cenário de trajetória bifurcada
    ├── Metricas_benchmark_dissertacao.rds                            # Lista consolidada com todas as tabelas de avaliação
    ├── Comparacao_final_metricas.csv                                 # Tabela final do benchmark com as configurações finais (PCA vs UMAP vs PHATE)
    └── Escalabilidade_metricas.csv                                   # Tabela com tempos e eficiência por número de células

---

## 3. Requisitos de Ambiente e Dependências
Para reproduzir integralmente os scripts e simulações, recomenda-se o seguinte ambiente computacional:

### Linguagens e Pacotes R (v >= 4.0)
- **Bioconductor:** `splatter`, `SingleCellExperiment`, `scran`, `scater`, `zellkonverter`, `Matrix`
- **Estatística e Métricas:** `mclust`, `stats`, `princurve`, `patchwork`
- **Ponte de Integração:** `reticulate`

### Linguagens e Pacotes Python (v >= 3.9)
- `scanpy` (modo *backed* para leitura de arquivos HDF5/AnnData)
- `scikit-learn` (`train_test_split`, `sklearn.metrics`, `sklearn.manifold`)
- `phate` (base C++/Python para integração com `phateR`)

---

## 4. Instruções de Execução
1. **Reprodutibilidade Estocástica:** Todos os scripts utilizam a semente aleatória global fixada em `set.seed(123)` (e equivalente em Python) para garantir a reprodutibilidade exata dos resultados.
2. **Ordem de Execução:**
  - Execute primeiramente o script Python `Amostragem_estratificada.py` para gerar o arquivo base.
  - Siga com as simulações executando sequencialmente os scripts iniciados por `Simulação_...` no R.
  - Para gerar os espaços latentes e métricas, execute os scripts de `Pré_processamento_embeddings_...`.
  - Por fim, consolide as visualizações com os scripts de `Métricas_...` e `Gráficos_tendência_...`.

---

## 5. Financiamento e Agradecimentos
Este trabalho foi realizado com o apoio de:
- **CAPES:** Código de Financiamento 001.
- **FAPESP:** Centro de Computação em Engenharia e Ciências (CCES - Proc. nº 2013/08293-7) e Cluster de Alto Desempenho COARACI (Proc. nº 2019/17874-0).

---

## 6. Citação
Se utilizar estes dados ou scripts em sua pesquisa, por favor cite:

ZELA, Alexis Marichi. **Métodos Estatísticos de Redução de Dimensionalidade em scRNA-seq: Fundamentos, Avaliação e Recomendações Práticas**. 2026. Dissertação (Mestrado em Estatística) - Instituto de Matemática, Estatística e Computação Científica (IMECC), Universidade Estadual de Campinas (UNICAMP), Campinas, 2026.
