#!/usr/bin/env Rscript

source(file.path("scripts", "load-fromage.R"))
load_fromage()

dir.create("results", showWarnings = FALSE)

lags <- c(1, 5, 10)

h0 <- c(ar1 = 0.4, ma1 = 0.3)
h1 <- c(ar1 = 0.6, ma1 = 0.5)
h2 <- c(ar1 = 0.2, ma1 = 0.1)

modelos_h1 <- list(h0, h1, h2)
h1_labels <- c(
  "ARMA(0.4,0.3)",
  "ARMA(0.6,0.5)",
  "ARMA(0.2,0.1)"
)

set.seed(2026)

resultados <- online_power_table_mm(
  lags = lags,
  h0_model = h0,
  h1_models = modelos_h1,
  n = 40,
  n_mc = 5,
  n_boot = 10,
  lag_q = 5,
  cores = 1,
  progress = FALSE,
  burn_in = 50
)

payload <- list(
  resultados = resultados,
  lags = lags,
  h0_model = h0,
  h1_models = modelos_h1,
  h1_labels = h1_labels,
  experiment = list(
    n = 40,
    n_mc = 5,
    n_boot = 10,
    alpha = 0.05,
    lag_q = 5,
    max.p = 2,
    max.q = 2,
    criterio = "aicc",
    order = c(1, 0, 1),
    burn_in = 50,
    seed = 2026
  )
)

saveRDS(payload, file = file.path("results", "small-experiment.rds"))

message("Saved results/small-experiment.rds")
