library(splatter)

sim <- readRDS("splatter_sim_10k.rds")

# Extrair o objeto de parâmetros original
params <- metadata(sim)$Params

# Mostrar os valores estimados
print(params)

# Salvar os parâmetros separadamente para poder modificá-los facilmente depois
saveRDS(params, "parametros_estimados_base.rds")

# Resultado

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
# 5  0.2, 0.2, 0.2, 0.2,...
# 
# Diff expr:
#   [PROBABILITY]             [DOWN PROB]              [LOCATION]
# 0.1, 0.1, 0.1, 0.1,...  0.5, 0.5, 0.5, 0.5,...  0.1, 0.1, 0.1, 0.1,...
# [SCALE]
# 0.4, 0.4, 0.4, 0.4,...
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
#   [FROM]                 [STEPS]                  [SKEW]
# 0, 0, 0, 0,...  100, 100, 100, 100,...  0.5, 0.5, 0.5, 0.5,...
# [Non-linear]          [Sigma Factor]
# 0.1                     0.8

# Mean, Library size, Exprs outliers, BCV son exactamente los valores puros
# estimados de los datos reales
