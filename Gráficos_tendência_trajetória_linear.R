#######################################################################
# Consolidação e Visualização de Resultados do Benchmark HVGs         #
#######################################################################

library(dplyr)
library(ggplot2)
library(tidyr)

# Identificar arquivos CSV
archivos_csv <- list.files(pattern = "^metricas_HVG_.*\\.csv$")

if (length(archivos_csv) == 0) {
  stop("Erro! Nenhum arquivo 'metricas_HVG_*.csv' foi encontrado. Verifique se os jobs do SLURM foram concluídos.")
}

message(paste("Foram encontrados", length(archivos_csv), "arquivos de métricas. Consolidando..."))

# Unir e ordenar dados
tabla_maestra <- do.call(rbind, lapply(archivos_csv, read.csv)) %>%
  arrange(HVG_Count, Method)

# Salvar Tabela Mestra CSV
write.csv(tabla_maestra, "TABLA_MAESTRA_EVALUACION_HVGS.csv", row.names = FALSE)
print(tabla_maestra)
message("\nTabela Mestra salva como 'TABLA_MAESTRA_EVALUACION_HVGS.csv'!")

# Gerar Gráficos de Tendência (Benchmark Final)
message("Gerando gráficos de tendência do Benchmark...")

# Transformar dados para o formato longo para o ggplot
tabla_larga <- tabla_maestra %>%
  pivot_longer(
    cols = c("Pseudotime_Spearman", "Global_Distance_Spearman", "Trustworthiness"),
    names_to = "Metrica",
    values_to = "Valor"
  )

p_tendencias <- ggplot(tabla_larga, aes(x = HVG_Count, y = Valor, color = Method, group = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~ Metrica, scales = "free_y") +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Impacto do Número de HVGs no Desempenho de DR",
    x = "Número de HVGs Selecionados",
    y = "Pontuação da Métrica",
    color = "Método DR"
  ) +
  theme_bw() +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

# Salvar figura final
ggsave("GRAFICO_BENCHMARK_HVG_TENDENCIAS.pdf", plot = p_tendencias, width = 12, height = 4.5)
message("Gráfico resumo salvo como 'GRAFICO_BENCHMARK_HVG_TENDENCIAS.pdf'!")