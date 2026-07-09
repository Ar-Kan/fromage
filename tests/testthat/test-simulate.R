test_that("conditional ARMA simulation returns the requested length", {
  set.seed(1)
  coefs <- c(ar1 = 0.4, ma1 = 0.3)

  simulated <- arima_structbreak_sim(
    n = 5,
    coefs = coefs,
    ordem = c(1, 0, 1),
    y = rnorm(20),
    eps = rnorm(20)
  )

  expect_type(simulated, "double")
  expect_length(simulated, 5)
  expect_false(anyNA(simulated))
})

test_that("structural break simulation returns both phases and innovations", {
  set.seed(2)
  coefs <- c(ar1 = 0.2, ma1 = 0.1)

  processo <- simula_processo(
    n = 20,
    fase1_coefs = coefs,
    fase2_coefs = coefs,
    order = c(1, 0, 1),
    burn_in = 20
  )

  expect_length(processo$fase1, 20)
  expect_length(processo$fase2, 20)
  expect_length(processo$eps_fase1, 20)
  expect_length(processo$eps_fase2, 20)
  expect_s3_class(processo$fit, "forecast_ARIMA")
})
