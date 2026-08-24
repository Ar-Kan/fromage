test_that("monitoring uses the documented coupled window", {
  set.seed(30)
  reference <- as.numeric(stats::arima.sim(list(ar = 0.3), n = 30))
  fit <- fit_qmm(reference, arma_candidates(1, 0))
  calibration <- calibrate_qmm(
    fit,
    horizons = 1:3,
    bootstrap_replicates = 4,
    portmanteau_lag = 4,
    seed = 31
  )
  new_data <- c(0.2, -0.1)
  monitored <- monitor_qmm(calibration, new_data)

  first_window <- c(reference[2:30], new_data[1])
  expected_first <- qmm_statistic(fit, first_window, portmanteau_lag = 4)

  expect_s3_class(monitored, "qmm_monitor")
  expect_equal(monitored$results$statistic[1], expected_first)
  expect_identical(
    monitored$results$signal,
    monitored$results$statistic > monitored$results$control_limit
  )
  expect_equal(signals(monitored), monitored$results[monitored$results$signal, ])
})

test_that("a monitoring object can accept another calibrated observation", {
  set.seed(32)
  reference <- stats::rnorm(30)
  fit <- fit_qmm(reference, arma_candidates(1, 0))
  calibration <- calibrate_qmm(
    fit,
    horizons = 1:3,
    bootstrap_replicates = 3,
    portmanteau_lag = 4,
    seed = 33
  )
  monitored <- monitor_qmm(calibration, c(0.1, 0.2))
  updated <- update(monitored, new_data = 0.3)

  expect_equal(updated$new_data, c(0.1, 0.2, 0.3))
  expect_equal(nrow(updated$results), 3)
  expect_equal(updated$series, c(reference, 0.1, 0.2, 0.3))
})

test_that("qmm_chart runs the complete single-series workflow", {
  simulated <- simulate_arma_change(
    reference_length = 30,
    monitoring_length = 3,
    reference_coefficients = c(ar1 = 0.3),
    monitoring_coefficients = c(ar1 = 0.6),
    order = c(1, 0),
    burn_in = 50,
    seed = 34
  )
  chart <- qmm_chart(
    simulated$series,
    reference_length = 30,
    candidates = arma_candidates(1, 0),
    bootstrap_replicates = 3,
    portmanteau_lag = 4,
    seed = 35
  )

  expect_s3_class(chart, "qmm_monitor")
  expect_equal(nrow(summary(chart)), 3)
  printed <- paste(capture.output(print(chart)), collapse = "\n")
  expect_match(printed, "Reference observations \\(Phase I\\): 30")
  expect_match(printed, "Monitoring observations \\(Phase II\\): 3")

  statistic_plot <- plot(chart)
  series_plot <- plot(chart, type = "series")
  overview_plot <- plot(chart, type = "overview")
  expect_s3_class(statistic_plot, "ggplot")
  expect_s3_class(series_plot, "ggplot")
  expect_s3_class(overview_plot, "ggplot")
  expect_equal(statistic_plot$labels$title, "Multimodel Q (QMM) control chart")
  expect_equal(series_plot$labels$subtitle, "The dashed line separates Phase I and Phase II")
  expect_equal(overview_plot$labels$title, "QMM Phase I and Phase II overview")
  expect_setequal(
    as.character(unique(overview_plot$data$panel)),
    c("Observed series", "QMM statistic")
  )
  expect_equal(
    sum(overview_plot$data$signal),
    2 * sum(chart$results$signal)
  )
})
