#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

input_file <- if (length(args) >= 1) args[[1]] else file.path("results", "small-experiment.rds")
output_file <- if (length(args) >= 2) args[[2]] else file.path("results", "power-plot.png")

suppressPackageStartupMessages({
  library(dplyr)
})

source(file.path("scripts", "load-fromage.R"))
load_fromage()

payload <- readRDS(input_file)
experiment <- payload$experiment

df_q <- cria_df_poder(
  resultados = payload$resultados,
  lags = payload$lags,
  h1_labels = payload$h1_labels,
  tipo_teste = "Q Multi-Modelo"
)

df_res <- cria_df_poder(
  resultados = payload$resultados,
  lags = payload$lags,
  h1_labels = payload$h1_labels,
  tipo_teste = "Res\u00edduo Auto-Arima"
)

df_plot <- dplyr::bind_rows(df_q, df_res)

info_experimento <- paste(
  sprintf("Modelo Fase 1 = ARMA(%s,%s)", payload$h0_model[["ar1"]], payload$h0_model[["ma1"]]),
  sprintf("Ordem = (%s,%s,%s)", experiment$order[1], experiment$order[2], experiment$order[3]),
  sprintf("MC = %s", experiment$n_mc),
  sprintf("Boot = %s", experiment$n_boot),
  sprintf("Alpha = %s", experiment$alpha),
  sprintf("Lag Q = %s", experiment$lag_q),
  sprintf("Lags = {%s}", paste(payload$lags, collapse = ", ")),
  sprintf("Modelos Candidatos Q: ARMA(max.p=%s,max.q=%s)", experiment$max.p, experiment$max.q),
  sprintf("Seed = %s", experiment$seed),
  sep = " | "
)

plot <- plot_power_curve(
  df_plot = df_plot,
  info_experimento = info_experimento,
  h1_labels = payload$h1_labels,
  alpha_line = experiment$alpha
)

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

ggplot2::ggsave(
  filename = output_file,
  plot = plot,
  device = "png",
  dpi = 600,
  width = 10,
  height = 6,
  units = "in",
  bg = "white"
)

message("Saved ", output_file)
