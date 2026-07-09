#!/usr/bin/env Rscript

source(file.path("scripts", "load-fromage.R"))
load_fromage()

dir.create("results", showWarnings = FALSE)

lags <- c(1, 5, 15, 25, 50, 99)

h0 <- c(ar1 = 0.4, ma1 = 0.3)
h1 <- c(ar1 = 0.6, ma1 = 0.5)
h2 <- c(ar1 = 0.8, ma1 = 0.7)
h3 <- c(ar1 = 0.2, ma1 = 0.1)
h4 <- c(ar1 = -0.4, ma1 = -0.3)

modelos_h1 <- list(h0, h1, h2, h3, h4)
h1_labels <- c(
  "ARMA(0.4,0.3)",
  "ARMA(0.6,0.5)",
  "ARMA(0.8,0.7)",
  "ARMA(0.2,0.1)",
  "ARMA(-0.4,-0.3)"
)

set.seed(2026)

resultados <- online_power_table_mm(
  lags = lags,
  h0_model = h0,
  h1_models = modelos_h1,
  n = 100,
  n_mc = 1000,
  n_boot = 500,
  alpha = 0.05,
  lag_q = 10,
  max.p = 2,
  max.q = 2,
  criterio = "aicc",
  order = c(1, 0, 1),
  burn_in = 1000,
  cores = 16,
  progress = TRUE
)

payload <- list(
  resultados = resultados,
  lags = lags,
  h0_model = h0,
  h1_models = modelos_h1,
  h1_labels = h1_labels,
  experiment = list(
    n = 100,
    n_mc = 1000,
    n_boot = 500,
    alpha = 0.05,
    lag_q = 10,
    max.p = 2,
    max.q = 2,
    criterio = "aicc",
    order = c(1, 0, 1),
    burn_in = 1000,
    cores = 16,
    seed = 2026
  )
)

saveRDS(payload, file = file.path("results", "paper-experiment.rds"))

message("Saved results/paper-experiment.rds")
