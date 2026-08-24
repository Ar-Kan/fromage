#' Fit the reference sample for a multimodel Q chart
#'
#' Fits each candidate `ARMA(p, q)` model to observations collected while the
#' process is believed to be in control. Models are compared using an
#' information criterion, and normalized relative-likelihood weights are
#' calculated within the candidate set.
#'
#' @param reference Numeric reference sample.
#' @param candidates Candidate models created by [arma_candidates()] or an
#'   equivalent data frame with columns `p` and `q`.
#' @param criterion Information criterion used for model weights: `"AICc"`,
#'   `"AIC"`, or `"BIC"`.
#' @param delta_max Largest information-criterion difference retained in the
#'   weighted set. The default, `Inf`, retains every successfully fitted model.
#' @param fit_method Method passed to [stats::arima()]. `"CSS-ML"` preserves
#'   the fitting convention used by the original prototype.
#'
#' @return A list with class `qmm_fit`. It contains the Phase I `reference`, a
#'   candidate table with fit diagnostics, criteria, differences, weights, and
#'   inclusion indicators, the corresponding fitted `stats::arima` objects,
#'   fitting settings, and the matched call.
#' @export
#'
#' @examples
#' set.seed(10)
#' reference <- as.numeric(stats::arima.sim(list(ar = 0.4, ma = 0.2), n = 60))
#' fit <- fit_qmm(reference, arma_candidates(1, 1))
#' fit
fit_qmm <- function(reference,
                    candidates = arma_candidates(),
                    criterion = c("AICc", "AIC", "BIC"),
                    delta_max = Inf,
                    fit_method = c("CSS-ML", "ML")) {
  reference <- validate_numeric_series(reference, "reference", min_length = 8L)
  candidates <- validate_candidates(candidates)
  criterion <- match.arg(criterion)
  fit_method <- match.arg(fit_method)

  if (!is.numeric(delta_max) || length(delta_max) != 1L || is.na(delta_max) || delta_max < 0) {
    stop("delta_max must be a single non-negative number or Inf.", call. = FALSE)
  }

  model_count <- nrow(candidates)
  models <- vector("list", model_count)
  status <- rep("failed", model_count)
  message <- rep("", model_count)
  log_likelihood <- rep(NA_real_, model_count)
  parameter_count <- rep(NA_integer_, model_count)
  aic <- aicc <- bic <- rep(Inf, model_count)

  for (i in seq_len(model_count)) {
    p <- candidates$p[i]
    q <- candidates$q[i]

    attempt <- tryCatch(
      capture_warnings(stats::arima(
        x = reference,
        order = c(p, 0L, q),
        include.mean = FALSE,
        method = fit_method
      )),
      error = function(error) error
    )

    if (inherits(attempt, "error")) {
      message[i] <- conditionMessage(attempt)
      next
    }

    model <- attempt$value
    if (!is.null(model$code) && model$code != 0L) {
      message[i] <- paste("Optimizer convergence code", model$code)
      next
    }

    ll <- tryCatch(stats::logLik(model), error = function(error) NULL)
    if (is.null(ll) || !is.finite(as.numeric(ll))) {
      message[i] <- "A finite maximized log-likelihood was not available."
      next
    }

    k <- length(model$coef) + 1L
    n <- length(reference)
    model_aic <- -2 * as.numeric(ll) + 2 * k
    model_aicc <- if (n - k - 1L > 0L) {
      model_aic + (2 * k * (k + 1L)) / (n - k - 1L)
    } else {
      Inf
    }

    models[[i]] <- model
    status[i] <- "fitted"
    message[i] <- paste(attempt$warnings, collapse = " | ")
    log_likelihood[i] <- as.numeric(ll)
    parameter_count[i] <- k
    aic[i] <- model_aic
    aicc[i] <- model_aicc
    bic[i] <- -2 * as.numeric(ll) + log(n) * k
  }

  criterion_values <- switch(criterion, AIC = aic, AICc = aicc, BIC = bic)
  finite <- is.finite(criterion_values)
  if (!any(finite)) {
    details <- paste(
      sprintf("ARMA(%d,%d): %s", candidates$p, candidates$q, message),
      collapse = "\n"
    )
    stop("No candidate model produced a finite ", criterion, ".\n", details, call. = FALSE)
  }

  delta <- rep(Inf, model_count)
  delta[finite] <- criterion_values[finite] - min(criterion_values[finite])
  included <- finite & delta <= delta_max
  relative_likelihood <- exp(-delta[included] / 2)
  weights <- numeric(model_count)
  weights[included] <- relative_likelihood / sum(relative_likelihood)

  model_table <- data.frame(
    model = sprintf("ARMA(%d,%d)", candidates$p, candidates$q),
    p = candidates$p,
    q = candidates$q,
    status = status,
    message = message,
    log_likelihood = log_likelihood,
    parameters = parameter_count,
    AIC = aic,
    AICc = aicc,
    BIC = bic,
    delta = delta,
    weight = weights,
    included = included,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      reference = reference,
      candidates = model_table,
      models = models,
      criterion = criterion,
      delta_max = delta_max,
      fit_method = fit_method,
      call = match.call()
    ),
    class = c("qmm_fit", "fromage_chart")
  )
}
