test_that("power results convert to tidy plotting data", {
  resultados <- list(
    resultados_q = list(
      li = matrix(c(0.01, 0.02), nrow = 1),
      media = matrix(c(0.05, 0.10), nrow = 1),
      ls = matrix(c(0.12, 0.18), nrow = 1)
    ),
    resultados_res = list(
      li = matrix(c(0.02, 0.03), nrow = 1),
      media = matrix(c(0.06, 0.11), nrow = 1),
      ls = matrix(c(0.13, 0.19), nrow = 1)
    )
  )

  df <- cria_df_poder(
    resultados = resultados,
    lags = c(1, 5),
    h1_labels = "ARMA(0.4,0.3)",
    tipo_teste = "Q Multi-Modelo"
  )

  expect_named(df, c("modelo", "Lag", "Poder", "Tipo", "LI", "LS"))
  expect_equal(nrow(df), 2)
  expect_equal(df$Lag, c(1, 5))
})
