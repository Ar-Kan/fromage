test_that("candidate grids use English columns and handle white noise", {
  candidates <- arma_candidates(max_p = 1, max_q = 1)

  expect_named(candidates, c("p", "q"))
  expect_equal(nrow(candidates), 3)
  expect_false(any(candidates$p == 0L & candidates$q == 0L))

  with_white_noise <- arma_candidates(0, 0, include_white_noise = TRUE)
  expect_equal(with_white_noise, data.frame(p = 0L, q = 0L))
})

test_that("candidate validation rejects duplicate and fractional orders", {
  reference <- seq_len(20)

  expect_error(
    fit_qmm(reference, data.frame(p = c(1, 1), q = c(0, 0))),
    "duplicate"
  )
  expect_error(
    fit_qmm(reference, data.frame(p = 0.5, q = 0)),
    "whole numbers"
  )
})

test_that("reference input must be a vector or univariate time series", {
  expect_error(
    fit_qmm(matrix(seq_len(20), ncol = 2), arma_candidates(1, 0)),
    "numeric vector or a univariate time series"
  )
})

test_that("fit_qmm reports information criteria and normalized weights", {
  set.seed(11)
  reference <- as.numeric(stats::arima.sim(list(ar = 0.45, ma = 0.2), n = 70))
  fit <- fit_qmm(reference, arma_candidates(1, 1))
  table <- summary(fit)

  expect_s3_class(fit, "qmm_fit")
  expect_named(
    table,
    c(
      "model", "p", "q", "status", "message", "log_likelihood",
      "parameters", "AIC", "AICc", "BIC", "delta", "weight", "included"
    )
  )
  expect_true(all(table$status == "fitted"))
  expect_equal(sum(table$weight), 1, tolerance = 1e-12)

  expected <- exp(-table$delta / 2)
  expected <- expected / sum(expected)
  expect_equal(table$weight, expected, tolerance = 1e-12)

  expected_aicc <- table$AIC +
    2 * table$parameters * (table$parameters + 1) /
      (length(reference) - table$parameters - 1)
  expect_equal(table$AICc, expected_aicc, tolerance = 1e-10)
})

test_that("delta_max retains a conditional subset and renormalizes weights", {
  set.seed(12)
  reference <- as.numeric(stats::arima.sim(list(ar = 0.5), n = 60))
  fit <- fit_qmm(reference, arma_candidates(2, 1), delta_max = 0)

  expect_equal(sum(fit$candidates$included), 1)
  expect_equal(sum(fit$candidates$weight), 1)
})
