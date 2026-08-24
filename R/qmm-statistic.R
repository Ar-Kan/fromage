#' Filter a series with fixed Phase I ARMA coefficients
#'
#' Applies an ARMA model to a new series without re-estimating its coefficients.
#' The returned residuals therefore measure how well the Phase I coefficient
#' values describe the evaluated window. ARMA(0,0) residuals are the input
#' series itself.
#'
#' @param x Numeric series to filter.
#' @param model A fitted `stats::arima` model whose coefficients remain fixed.
#' @param p Autoregressive order of `model`.
#' @param q Moving-average order of `model`.
#'
#' @return A numeric residual vector with the same length as `x`.
#' @keywords internal
fixed_arma_residuals <- function(x, model, p, q) {
  if (p == 0L && q == 0L) {
    return(x)
  }

  fixed_fit <- stats::arima(
    x = x,
    order = c(p, 0L, q),
    include.mean = FALSE,
    fixed = model$coef,
    transform.pars = FALSE,
    method = "ML"
  )

  as.numeric(fixed_fit$residuals)
}

#' Calculate the model-specific components of the QMM statistic
#'
#' For every retained candidate, filters the evaluated series with fixed Phase
#' I coefficients and calculates a Ljung-Box statistic with `fitdf = 0`. The
#' candidate weights are copied from the Phase I fit; this helper does not
#' combine the components.
#'
#' @param fit A validated `qmm_fit` object.
#' @param x Numeric series or coupled monitoring window to evaluate.
#' @param portmanteau_lag Number of residual autocorrelations used in each
#'   Ljung-Box statistic.
#'
#' @return A data frame with one retained candidate per row and columns
#'   `model`, `statistic`, and `weight`.
#' @keywords internal
qmm_components <- function(fit, x, portmanteau_lag) {
  included <- which(fit$candidates$included)
  statistics <- numeric(length(included))

  for (j in seq_along(included)) {
    i <- included[j]
    residuals <- fixed_arma_residuals(
      x,
      fit$models[[i]],
      fit$candidates$p[i],
      fit$candidates$q[i]
    )
    if (anyNA(residuals) || any(!is.finite(residuals))) {
      stop("Residual filtering failed for ", fit$candidates$model[i], ".", call. = FALSE)
    }

    statistics[j] <- as.numeric(stats::Box.test(
      residuals,
      lag = portmanteau_lag,
      type = "Ljung-Box",
      fitdf = 0L
    )$statistic)
  }

  data.frame(
    model = fit$candidates$model[included],
    statistic = statistics,
    weight = fit$candidates$weight[included],
    stringsAsFactors = FALSE
  )
}

#' Calculate the multimodel Q statistic
#'
#' Calculates a Ljung-Box statistic for the residuals produced by each retained
#' candidate model and returns their information-criterion-weighted average.
#' Model coefficients and weights remain fixed at their reference-sample values.
#'
#' @param fit A fitted reference sample created by [fit_qmm()].
#' @param x Numeric series to evaluate. For monitoring, this is the coupled
#'   fixed-length window.
#' @param portmanteau_lag Number of residual autocorrelations included in each
#'   Ljung-Box statistic.
#'
#' @return A single numeric value.
#' @export
qmm_statistic <- function(fit, x, portmanteau_lag = 10L) {
  validate_qmm_fit(fit)
  x <- validate_numeric_series(x, "x", min_length = 3L)
  portmanteau_lag <- validate_whole_number(portmanteau_lag, "portmanteau_lag")
  if (portmanteau_lag >= length(x)) {
    stop("portmanteau_lag must be smaller than the evaluated series length.", call. = FALSE)
  }

  components <- qmm_components(fit, x, portmanteau_lag)
  as.numeric(sum(components$statistic * components$weight))
}
