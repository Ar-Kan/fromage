#' Monitor new observations with a calibrated multimodel Q chart
#'
#' Evaluates new observations in their arrival order. At horizon `L`, the
#' window contains the final `T - L` reference observations followed by the
#' first `L` new observations, where `T` is the reference-sample length.
#'
#' @param calibration Control limits created by [calibrate_qmm()].
#' @param new_data Numeric observations in monitoring order.
#'
#' @return A list with class `qmm_monitor` containing the calibration, Phase II
#'   observations, horizon-specific statistics and limits, signal indicators,
#'   the complete series, the Phase I length, and the matched call.
#' @export
monitor_qmm <- function(calibration, new_data) {
  validate_qmm_calibration(calibration)
  new_data <- validate_numeric_series(new_data, "new_data", min_length = 1L)
  fit <- calibration$fit
  reference_length <- length(fit$reference)

  if (length(new_data) > reference_length) {
    stop(
      "new_data cannot be longer than the reference sample under the QMM construction.",
      call. = FALSE
    )
  }

  horizons <- seq_along(new_data)
  limit_rows <- match(horizons, calibration$limits$horizon)
  if (anyNA(limit_rows)) {
    missing <- horizons[is.na(limit_rows)]
    stop(
      "The calibration does not contain horizon(s): ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  statistics <- numeric(length(new_data))
  for (horizon in horizons) {
    window <- if (horizon < reference_length) {
      c(
        fit$reference[(horizon + 1L):reference_length],
        new_data[seq_len(horizon)]
      )
    } else {
      new_data
    }
    statistics[horizon] <- qmm_statistic(
      fit,
      window,
      portmanteau_lag = calibration$portmanteau_lag
    )
  }

  control_limits <- calibration$limits$control_limit[limit_rows]
  results <- data.frame(
    horizon = horizons,
    observation = new_data,
    statistic = statistics,
    control_limit = control_limits,
    signal = statistics > control_limits
  )

  structure(
    list(
      calibration = calibration,
      new_data = new_data,
      results = results,
      series = c(fit$reference, new_data),
      reference_length = reference_length,
      call = match.call()
    ),
    class = c("qmm_monitor", "fromage_chart")
  )
}

#' Fit, calibrate, and monitor one temporal series
#'
#' A convenient complete workflow for a series containing a reference period
#' followed by new observations. Use [fit_qmm()], [calibrate_qmm()], and
#' [monitor_qmm()] separately when each intermediate result should be inspected.
#'
#' @param x Complete numeric series.
#' @param reference_length Number of initial observations treated as the
#'   in-control reference sample.
#' @inheritParams fit_qmm
#' @inheritParams calibrate_qmm
#'
#' @return A `qmm_monitor` object as described by [monitor_qmm()].
#' @export
#'
#' @examples
#' simulated <- simulate_arma_change(
#'   reference_length = 80,
#'   monitoring_length = 30,
#'   reference_coefficients = c(ar1 = 0.4, ma1 = 0.3),
#'   monitoring_coefficients = c(ar1 = 0.9, ma1 = 0.85),
#'   order = c(1, 1),
#'   seed = 2026
#' )
#'
#' chart <- qmm_chart(
#'   simulated$series,
#'   reference_length = 80,
#'   candidates = arma_candidates(1, 1),
#'   bootstrap_replicates = 30,
#'   portmanteau_lag = 10,
#'   seed = 99
#' )
#'
#' chart
#' signals(chart)
#' plot(chart, type = "overview")
qmm_chart <- function(x,
                      reference_length,
                      candidates = arma_candidates(),
                      criterion = c("AICc", "AIC", "BIC"),
                      delta_max = Inf,
                      fit_method = c("CSS-ML", "ML"),
                      bootstrap_replicates = 500L,
                      confidence_level = 0.95,
                      portmanteau_lag = 10L,
                      seed = NULL,
                      max_draw_attempts = 1000L,
                      keep_bootstrap = TRUE) {
  x <- validate_numeric_series(x, "x", min_length = 3L)
  reference_length <- validate_whole_number(reference_length, "reference_length", minimum = 2L)
  if (reference_length >= length(x)) {
    stop("reference_length must leave at least one new observation.", call. = FALSE)
  }

  reference <- x[seq_len(reference_length)]
  new_data <- x[(reference_length + 1L):length(x)]
  if (length(new_data) > reference_length) {
    stop(
      "The number of new observations cannot exceed reference_length under the QMM construction.",
      call. = FALSE
    )
  }

  fit <- fit_qmm(
    reference = reference,
    candidates = candidates,
    criterion = criterion,
    delta_max = delta_max,
    fit_method = fit_method
  )
  calibration <- calibrate_qmm(
    fit = fit,
    horizons = seq_along(new_data),
    bootstrap_replicates = bootstrap_replicates,
    confidence_level = confidence_level,
    portmanteau_lag = portmanteau_lag,
    seed = seed,
    max_draw_attempts = max_draw_attempts,
    keep_bootstrap = keep_bootstrap
  )
  monitored <- monitor_qmm(calibration, new_data)
  monitored$call <- match.call()
  monitored
}

#' Extract signalled observations
#'
#' @param object A monitoring result.
#' @param ... Additional arguments passed to an object-specific method.
#'
#' @return A data frame containing only signalled horizons.
#' @export
signals <- function(object, ...) {
  UseMethod("signals")
}

#' @rdname signals
#' @export
signals.qmm_monitor <- function(object, ...) {
  object$results[object$results$signal, , drop = FALSE]
}

#' Add Phase II observations to a QMM monitoring result
#'
#' Re-evaluates the accumulated Phase II observations using the limits already
#' stored in the calibration. It does not refit the Phase I models or rerun the
#' bootstrap. Every resulting horizon must therefore be present in the
#' calibration.
#'
#' @param object A `qmm_monitor` object.
#' @param new_data One or more additional Phase II observations, in arrival
#'   order.
#' @param ... Unused.
#'
#' @return A new `qmm_monitor` object containing the previous and additional
#'   Phase II observations.
#' @method update qmm_monitor
#' @export
update.qmm_monitor <- function(object, new_data, ...) {
  additional <- validate_numeric_series(new_data, "new_data", min_length = 1L)
  monitor_qmm(object$calibration, c(object$new_data, additional))
}
