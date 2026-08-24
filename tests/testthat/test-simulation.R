test_that("ARMA change simulation is reproducible and returns both periods", {
  first <- simulate_arma_change(
    reference_length = 20,
    monitoring_length = 5,
    reference_coefficients = c(ar1 = 0.4, ma1 = 0.2),
    monitoring_coefficients = c(ar1 = 0.6, ma1 = 0.3),
    order = c(1, 1),
    burn_in = 30,
    seed = 40
  )
  second <- simulate_arma_change(
    reference_length = 20,
    monitoring_length = 5,
    reference_coefficients = c(ar1 = 0.4, ma1 = 0.2),
    monitoring_coefficients = c(ar1 = 0.6, ma1 = 0.3),
    order = c(1, 1),
    burn_in = 30,
    seed = 40
  )

  expect_equal(first, second)
  expect_length(first$series, 25)
  expect_equal(first$series, c(first$reference, first$monitoring))
  expect_length(first$reference_innovations, 20)
  expect_length(first$monitoring_innovations, 5)
})

test_that("simulation rejects nonstationary and noninvertible coefficients", {
  expect_error(
    simulate_arma_change(20, 5, c(ar1 = 1.1), c(ar1 = 0.2), order = c(1, 0)),
    "stationary"
  )
  expect_error(
    simulate_arma_change(20, 5, c(ma1 = 1.1), c(ma1 = 0.2), order = c(0, 1)),
    "invertible"
  )
})

test_that("simulation uses zero pre-sample values when burn-in is zero", {
  simulated <- simulate_arma_change(
    reference_length = 3,
    monitoring_length = 1,
    reference_coefficients = c(ar1 = 0.4),
    monitoring_coefficients = c(ar1 = 0.4),
    order = c(1, 0),
    burn_in = 0,
    seed = 41
  )

  expect_equal(simulated$series[1], simulated$reference_innovations[1])
})
