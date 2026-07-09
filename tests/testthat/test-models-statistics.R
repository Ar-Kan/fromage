test_that("candidate grid excludes white noise by default", {
  grid <- arma_candidates(max.p = 1, max.q = 1)

  expect_named(grid, c("p", "q"))
  expect_false(any(grid$p == 0 & grid$q == 0))
  expect_equal(nrow(grid), 3)
})

test_that("model fitting and multimodel Q return usable outputs", {
  set.seed(3)
  serie <- as.numeric(stats::arima.sim(
    model = list(ar = 0.3, ma = 0.2),
    n = 50
  ))

  fit <- ajustar_modelos(
    serie = serie,
    modelos_df = arma_candidates(max.p = 1, max.q = 1),
    criterio = "aicc"
  )

  expect_named(fit, c("modelos", "pesos", "ics", "modelos_df"))
  expect_equal(sum(fit$pesos), 1, tolerance = 1e-8)

  q_stat <- multimodel_Q(
    serie = serie,
    resultados_mmi = fit,
    lag_q = 5
  )

  expect_type(q_stat, "double")
  expect_length(q_stat, 1)
  expect_true(is.finite(q_stat))
})
