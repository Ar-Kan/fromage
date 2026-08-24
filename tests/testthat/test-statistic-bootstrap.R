test_that("qmm_statistic is the weighted sum of candidate statistics", {
  set.seed(20)
  reference <- as.numeric(stats::arima.sim(list(ar = 0.35, ma = 0.2), n = 55))
  fit <- fit_qmm(reference, arma_candidates(1, 1))

  components <- fromage:::qmm_components(fit, reference, portmanteau_lag = 5)
  observed <- qmm_statistic(fit, reference, portmanteau_lag = 5)

  expect_equal(observed, sum(components$statistic * components$weight))
  expect_true(is.finite(observed))
})

test_that("portmanteau lag must fit inside the evaluated window", {
  set.seed(21)
  reference <- stats::rnorm(30)
  fit <- fit_qmm(reference, arma_candidates(1, 0))

  expect_error(qmm_statistic(fit, reference, 30), "smaller")
})

test_that("calibration is reproducible and restores the caller RNG state", {
  set.seed(22)
  reference <- as.numeric(stats::arima.sim(list(ar = 0.4), n = 35))
  fit <- fit_qmm(reference, arma_candidates(1, 0))

  set.seed(999)
  state_before <- .Random.seed
  first <- calibrate_qmm(
    fit,
    horizons = 1:3,
    bootstrap_replicates = 5,
    portmanteau_lag = 4,
    seed = 100
  )
  expect_identical(.Random.seed, state_before)

  second <- calibrate_qmm(
    fit,
    horizons = 1:3,
    bootstrap_replicates = 5,
    portmanteau_lag = 4,
    seed = 100
  )

  expect_equal(first$bootstrap, second$bootstrap)
  expect_equal(first$limits, second$limits)
  expect_true(all(first$diagnostics$rejected_draws >= 0))
  expect_true(all(first$selected_models %in% which(fit$candidates$included)))
})

test_that("calibration validates the supported horizon range", {
  set.seed(23)
  reference <- stats::rnorm(20)
  fit <- fit_qmm(reference, arma_candidates(1, 0))

  expect_error(
    calibrate_qmm(fit, 21, bootstrap_replicates = 2, portmanteau_lag = 3),
    "no greater"
  )
})

test_that("bootstrap sampling handles one retained model at a non-first row", {
  set.seed(24)
  reference <- stats::rnorm(30)
  fit <- fit_qmm(
    reference,
    data.frame(p = c(0L, 1L), q = c(0L, 0L))
  )
  fit$candidates$included <- c(FALSE, TRUE)
  fit$candidates$weight <- c(0, 1)

  calibration <- calibrate_qmm(
    fit,
    horizons = 1:2,
    bootstrap_replicates = 3,
    portmanteau_lag = 4,
    seed = 25
  )

  expect_true(all(calibration$selected_models == 2L))
})

test_that("white-noise candidates can be calibrated and monitored", {
  set.seed(26)
  reference <- stats::rnorm(25)
  fit <- fit_qmm(
    reference,
    arma_candidates(0, 0, include_white_noise = TRUE)
  )
  calibration <- calibrate_qmm(
    fit,
    horizons = 1:2,
    bootstrap_replicates = 3,
    portmanteau_lag = 3,
    seed = 27
  )
  monitored <- monitor_qmm(calibration, c(0.1, -0.1))

  expect_true(all(is.finite(monitored$results$statistic)))
  expect_true(all(is.finite(monitored$results$control_limit)))
})
