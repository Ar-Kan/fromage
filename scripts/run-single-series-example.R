#!/usr/bin/env Rscript
# saves a reproducible result and overview plot under `results/`

source(file.path("scripts", "load-fromage.R"))
load_fromage()

simulated <- simulate_arma_change(
  reference_length = 100,
  monitoring_length = 40,
  reference_coefficients = c(ar1 = 0.4, ma1 = 0.3),
  monitoring_coefficients = c(ar1 = 0.8, ma1 = 0.7),
  order = c(1, 1),
  seed = 2026
)

chart <- qmm_chart(
  x = simulated$series,
  reference_length = 100,
  candidates = arma_candidates(max_p = 1, max_q = 1),
  bootstrap_replicates = 199,
  portmanteau_lag = 10,
  seed = 99
)

dir.create("results", showWarnings = FALSE)
saveRDS(chart, file.path("results", "single-series-example.rds"))
ggplot2::ggsave(
  filename = file.path("results", "single-series-example.png"),
  plot = plot(chart, type = "overview"),
  width = 7,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)

print(summary(chart))
print(signals(chart))
message("Saved the example result and chart in results/.")
