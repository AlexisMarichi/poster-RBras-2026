#########################################
# Simulação de Trajetória Sequencial    #
#########################################

set.seed(123) # Mantemos a mesma semente para comparabilidade

library(splatter)
library(SingleCellExperiment)

message("Carregando o molde de parâmetros base...")
# params <- readRDS("parametros_estimados_base.rds")
params_base <- readRDS("parametros_estimados_base.rds")

# CONFIGURAÇÃO DA TRAJETÓRIA LINEAR 

# Todos estes parâmetros abaixo devem ser alterados antes de realizar a simulação 

# Criar um objeto completamente novo (por padrão nasce com 1 grupo)
params <- newSplatParams()

# Copiar os parâmetros estimados globais e aplicar os da trajetória
params <- setParams(
  params,
  # Copiar as estimativas globais do arquivo RDS 
  nGenes = getParam(params_base, "nGenes"),
  batchCells  = getParam(params_base, "batchCells"),
  seed = getParam(params_base, "seed"), 
  mean.rate = getParam(params_base, "mean.rate"),
  mean.shape = getParam(params_base, "mean.shape"),
  lib.loc = getParam(params_base, "lib.loc"),
  lib.scale = getParam(params_base, "lib.scale"),
  out.prob = getParam(params_base, "out.prob"),
  out.facLoc = getParam(params_base, "out.facLoc"),
  out.facScale = getParam(params_base, "out.facScale"),
  bcv.common = getParam(params_base, "bcv.common"),
  bcv.df = getParam(params_base, "bcv.df"),
  dropout.mid = getParam(params_base, "dropout.mid"),
  dropout.shape = getParam(params_base, "dropout.shape"),
  
  # CONFIGURAÇÃO DE TRAJETÓRIA SEQUENCIAL
  
  # Dividimos as células equitativamente em 4 segmentos temporais
  group.prob = c(0.25, 0.25, 0.25, 0.25),
  
  # Nós de origem: 
  # O Path 1 vem da origem (0)
  # O Path 2 nasce ao final do Path 1
  # O Path 3 nasce ao final do Path 2... etc.
  path.from = c(0, 1, 2, 3),
  # INTENSIDADE E RUÍDO DA TRAJETÓRIA
  de.prob    = 0.5,   # Valor padrão 0.1; High DE: 50% dos genes marcam a mudança de estado
  de.facLoc  = 0.5,   # Magnitude da mudança, valor padrão 0.1
  path.sigmaFac = 0.1, # Ruído baixo para manter a curva definida, valor padrão 0.8
  # path.nonlinearProb = 0 # valor padrão 0.1 e elimina distorções não lineares para um controle estrito do benchmark
)

# params <- setParam(params, "group.prob", 1)
# 
# PARÂMETROS ESPECÍFICOS DE PATHS
# params <- setParam(params, "path.from", 0)             # NÓ DE ORIGEM: 0 indica que inicia a partir da origem, por padrão
# params <- setParam(params, "path.nSteps", 100)         # PASSOS: Resolução da trajetória (100 passos discretos), por padrão
# params <- setParam(params, "path.skew", 0.5)           # ASSIMETRIA (SKEW): 0.5 = distribuição uniforme de células no caminho, por padrão
# params <- setParam(params, "path.nonlinearProb", 0.1) # NÃO LINEARIDADE: Grau de curva na resposta gênica, por padrão
# params <- setParam(params, "path.sigmaFac", 0.8)    # RUÍDO: Variabilidade aleatória sobre a tendência do gene, por padrão
# 
# AJUSTANDO OUTROS PARÂMETROS PARA QUE OS PATHS SEJAM BEM SIMULADOS
# 
# params <- setParam(params, "de.prob", 0.1)        # 10% dos genes mudam ao longo do caminho
# params <- setParam(params, "de.downProb", 0.5)    # 50% de probabilidade de subexpressão
# params <- setParam(params, "de.facLoc", 0.5)      # Média do fator DE, padrão 0.1
# params <- setParam(params, "de.facScale", 0.4)    # Desvio padrão do fator DE

# Verificação no console para confirmar a alteração
message("\n--- Verificando alterações no objeto params ---")
print(params)

# A Params object of class SplatParams
# Parameters can be (estimable) or [not estimable], 'Default' or  'NOT DEFAULT'
# Secondary parameters are usually set during simulation
# 
# Global:
#   (GENES)  (CELLS)   [SEED]
# 42206    10000   722449
# 
# 29 additional parameters
# 
# Batches:
#   [BATCHES]  [BATCH CELLS]     [Location]        [Scale]       [Remove]
# 1          10000            0.1            0.1          FALSE
# 
# Mean:
#   (RATE)            (SHAPE)
# 1.07655534549026  0.222735704297571
# 
# Library size:
#   (LOCATION)            (SCALE)             (Norm)
# 9.75586695405308  0.838332562986004              FALSE
# 
# Exprs outliers:
#   (PROBABILITY)            (LOCATION)               (SCALE)
# 0.000529033883360625      8.40837810203532     0.627160811415717
# 
# Groups:
#   [GROUPS]           [GROUP PROBS]
# 4  0.25, 0.25, 0.25, 0.25
# 
# Diff expr:
#   [PROBABILITY]    [Down Prob]     [LOCATION]        [Scale]
# 0.5            0.5            0.5            0.4
# 
# BCV:
#   (COMMON DISP)              (DOF)
# 0.265781616479207   8.26184462563436
# 
# Dropout:
#   [Type]           (MIDPOINT)              (SHAPE)
# none  -0.0999198989660964    -1.02193333868232
# 
# Paths:
#   [FROM]         [Steps]          [Skew]    [Non-linear]  [SIGMA FACTOR]
# 0, 1, 2, 3             100             0.5             0.1             0.1

# SIMULAÇÃO COM TRAJETÓRIA LINEAR

message("\nSimulando dados com trajetória linear (method = 'paths')...")
set.seed(123)

sim_path <- splatSimulate(
  params,
  nGenes = 10000,
  batchCells = 1000,
  method = "paths",
  verbose = TRUE
)

# Salvamos o arquivo resultante
nombre_salida <- paste0("splatter_sim_10k_linear_path", ".rds")
saveRDS(sim_path, nombre_salida)

message(paste("Processo concluído! Arquivo salvo como:", nombre_salida))

params_simulacion <- metadata(sim_path)$Params
message("\n=== PARÂMETROS USADOS NA SIMULAÇÃO FINAL ===")
print(params_simulacion)

# A Params object of class SplatParams
# Parameters can be (estimable) or [not estimable], 'Default' or  'NOT DEFAULT'
# Secondary parameters are usually set during simulation
# 
# Global:
#   (Genes)  (CELLS)   [SEED]
# 10000     1000   722449
# 
# 29 additional parameters
# 
# Batches:
#   [BATCHES]  [BATCH CELLS]     [Location]        [Scale]       [Remove]
# 1           1000            0.1            0.1          FALSE
# 
# Mean:
#   (RATE)            (SHAPE)
# 1.07655534549026  0.222735704297571
# 
# Library size:
#   (LOCATION)            (SCALE)             (Norm)
# 9.75586695405308  0.838332562986004              FALSE
# 
# Exprs outliers:
#   (PROBABILITY)            (LOCATION)               (SCALE)
# 0.000529033883360625      8.40837810203532     0.627160811415717
# 
# Groups:
#   [GROUPS]           [GROUP PROBS]
# 4  0.25, 0.25, 0.25, 0.25
# 
# Diff expr:
#   [PROBABILITY]         [DOWN PROB]          [LOCATION]             [SCALE]
# 0.5, 0.5, 0.5, 0.5  0.5, 0.5, 0.5, 0.5  0.5, 0.5, 0.5, 0.5  0.4, 0.4, 0.4, 0.4
# 
# BCV:
#   (COMMON DISP)              (DOF)
# 0.265781616479207   8.26184462563436
# 
# Dropout:
#   [Type]           (MIDPOINT)              (SHAPE)
# none  -0.0999198989660964    -1.02193333868232
# 
# Paths:
#   [FROM]             [STEPS]              [SKEW]        [Non-linear]
# 0, 1, 2, 3  100, 100, 100, 100  0.5, 0.5, 0.5, 0.5                 0.1
# [SIGMA FACTOR]
# 0.1