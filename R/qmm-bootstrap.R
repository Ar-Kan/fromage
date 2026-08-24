#' Draw admissible ARMA coefficients for one bootstrap replicate
#'
#' Draws a coefficient vector from the multivariate normal approximation
#' defined by a fitted ARMA model's coefficient estimates and covariance
#' matrix. Draws that do not define a stationary and invertible model are
#' rejected. This helper is called by [bootstrap_qmm_once()].
#'
#' @param model A fitted `stats::arima` object from the Phase I sample.
#' @param p Autoregressive order of `model`.
#' @param q Moving-average order of `model`.
#' @param max_attempts Maximum number of coefficient draws before reporting an
#'   error.
#'
#' @return A list with `coefficients`, the accepted named numeric vector, and
#'   `attempts`, the number of draws required. An ARMA(0,0) model returns an
#'   empty coefficient vector on the first attempt.
#' @keywords internal
draw_arma_coefficients <- function(model, p, q, max_attempts) {
  if (length(model$coef) == 0L) {
    return(list(coefficients = numeric(), attempts = 1L))
  }
  if (is.null(model$var.coef) || anyNA(model$var.coef) || any(!is.finite(model$var.coef))) {
    stop("The fitted coefficient covariance matrix is not finite.", call. = FALSE)
  }

  for (attempt in seq_len(max_attempts)) {
    draw <- as.numeric(mvtnorm::rmvnorm(1L, mean = model$coef, sigma = model$var.coef))
    names(draw) <- names(model$coef)
    if (arma_is_admissible(draw, p, q)) {
      return(list(coefficients = draw, attempts = attempt))
    }
  }

  stop(
    "Unable to draw stationary and invertible ARMA coefficients after ",
    max_attempts,
    " attempts.",
    call. = FALSE
  )
}

#' Simulate an ARMA continuation conditional on existing histories
#'
#' Continues an ARMA recursion using fixed coefficients. The most recent
#' observed values supply the autoregressive history, and the supplied
#' innovations supply the moving-average history. New Gaussian innovations
#' are generated with standard deviation `innovation_sd`. This is the
#' conditional simulation step used by the parametric bootstrap.
#'
#' @param n Number of new observations to simulate.
#' @param coefficients Named ARMA coefficient vector.
#' @param p Autoregressive order.
#' @param q Moving-average order.
#' @param observations Numeric history observed before the continuation.
#' @param innovations Numeric innovation history before the continuation. In
#'   the feasible package workflow, this is the selected model's estimated
#'   residual history.
#' @param innovation_sd Standard deviation of the new Gaussian innovations.
#'
#' @return A numeric vector containing `n` simulated observations. The supplied
#'   histories are used as initial conditions and are not included in the
#'   returned vector.
#' @keywords internal
simulate_arma_continuation <- function(n,
                                       coefficients,
                                       p,
                                       q,
                                       observations,
                                       innovations,
                                       innovation_sd) {
  ar <- if (p > 0L) coefficients[paste0("ar", seq_len(p))] else numeric()
  ma <- if (q > 0L) coefficients[paste0("ma", seq_len(q))] else numeric()
  new_innovations <- stats::rnorm(n, sd = innovation_sd)
  new_observations <- numeric(n)
  observation_history <- as.numeric(observations)
  innovation_history <- as.numeric(innovations)

  for (time in seq_len(n)) {
    ar_part <- if (p > 0L) {
      sum(ar * rev(utils::tail(observation_history, p)))
    } else {
      0
    }
    ma_part <- if (q > 0L) {
      sum(ma * rev(utils::tail(innovation_history, q)))
    } else {
      0
    }

    new_observations[time] <- ar_part + ma_part + new_innovations[time]
    observation_history <- c(observation_history, new_observations[time])
    innovation_history <- c(innovation_history, new_innovations[time])
  }

  new_observations
}

#' Generate one QMM bootstrap statistic
#'
#' Samples one retained candidate according to its Phase I weight, draws
#' admissible coefficients, simulates an in-control continuation, constructs
#' the coupled window for the requested horizon, and evaluates its QMM
#' statistic. This is one replicate of the procedure orchestrated by
#' [calibrate_qmm()].
#'
#' @param fit A validated `qmm_fit` object.
#' @param horizon Monitoring horizon for the coupled window.
#' @param portmanteau_lag Number of residual autocorrelations used in the
#'   Ljung-Box statistics.
#' @param max_draw_attempts Maximum coefficient draws allowed when searching
#'   for a stationary and invertible draw.
#'
#' @return A named numeric vector with `statistic`, `draw_attempts`, and
#'   `selected_model`. The last value is the row number of the sampled model in
#'   `fit$candidates`.
#' @keywords internal
bootstrap_qmm_once <- function(fit, horizon, portmanteau_lag, max_draw_attempts) {
  included <- which(fit$candidates$included)
  selected_position <- sample.int(
    length(included),
    size = 1L,
    prob = fit$candidates$weight[included]
  )
  selected <- included[selected_position]
  model <- fit$models[[selected]]
  p <- fit$candidates$p[selected]
  q <- fit$candidates$q[selected]
  draw <- draw_arma_coefficients(model, p, q, max_draw_attempts)

  residual_history <- as.numeric(model$residuals)
  residual_history <- residual_history[is.finite(residual_history)]
  if (length(residual_history) < q) {
    stop("The selected model does not provide enough finite residuals.", call. = FALSE)
  }

  simulated <- simulate_arma_continuation(
    n = horizon,
    coefficients = draw$coefficients,
    p = p,
    q = q,
    observations = fit$reference,
    innovations = residual_history,
    innovation_sd = sqrt(model$sigma2)
  )

  reference_length <- length(fit$reference)
  window <- if (horizon < reference_length) {
    c(fit$reference[(horizon + 1L):reference_length], simulated)
  } else {
    simulated
  }

  c(
    statistic = qmm_statistic(fit, window, portmanteau_lag),
    draw_attempts = draw$attempts,
    selected_model = selected
  )
}

#' Calibrate multimodel Q control limits
#'
#' Uses a parametric bootstrap to estimate an upper control limit for every
#' requested monitoring horizon. Each bootstrap replicate samples a candidate
#' model using its reference-sample weight, samples admissible coefficients,
#' simulates an in-control continuation, and evaluates the coupled window.
#'
#' @param fit A reference-sample fit created by [fit_qmm()].
#' @param horizons Whole-number monitoring horizons. Every value must be no
#'   greater than the reference-sample length.
#' @param bootstrap_replicates Number of bootstrap replicates per horizon.
#' @param confidence_level Pointwise control-limit probability, such as `0.95`.
#' @param portmanteau_lag Number of residual autocorrelations used by
#'   [qmm_statistic()].
#' @param seed Optional random-number seed. The caller's random-number state is
#'   restored when the function returns.
#' @param max_draw_attempts Maximum coefficient draws allowed per replicate.
#' @param keep_bootstrap Whether to retain the bootstrap statistics in the
#'   returned object.
#'
#' @return A list with class `qmm_calibration`. It contains the Phase I `fit`,
#'   horizon-specific `limits`, optional bootstrap draws and selected-model
#'   indices, coefficient-draw diagnostics, calibration settings, and the
#'   matched call.
#' @export
calibrate_qmm <- function(fit,
                          horizons,
                          bootstrap_replicates = 500L,
                          confidence_level = 0.95,
                          portmanteau_lag = 10L,
                          seed = NULL,
                          max_draw_attempts = 1000L,
                          keep_bootstrap = TRUE) {
  validate_qmm_fit(fit)
  bootstrap_replicates <- validate_whole_number(
    bootstrap_replicates,
    "bootstrap_replicates"
  )
  confidence_level <- validate_probability(confidence_level, "confidence_level")
  portmanteau_lag <- validate_whole_number(portmanteau_lag, "portmanteau_lag")
  max_draw_attempts <- validate_whole_number(max_draw_attempts, "max_draw_attempts")

  if (!is.numeric(horizons) || length(horizons) == 0L || anyNA(horizons) ||
      any(!is.finite(horizons)) || any(horizons != floor(horizons)) || any(horizons < 1L)) {
    stop("horizons must contain positive whole numbers.", call. = FALSE)
  }
  horizons <- sort(unique(as.integer(horizons)))
  reference_length <- length(fit$reference)
  if (any(horizons > reference_length)) {
    stop("Every horizon must be no greater than the reference-sample length.", call. = FALSE)
  }
  if (portmanteau_lag >= reference_length) {
    stop("portmanteau_lag must be smaller than the reference-sample length.", call. = FALSE)
  }

  bootstrap <- matrix(
    NA_real_,
    nrow = bootstrap_replicates,
    ncol = length(horizons),
    dimnames = list(NULL, paste0("horizon_", horizons))
  )
  attempts <- matrix(NA_integer_, nrow = bootstrap_replicates, ncol = length(horizons))
  selected_models <- matrix(NA_integer_, nrow = bootstrap_replicates, ncol = length(horizons))

  with_seed(seed, {
    for (column in seq_along(horizons)) {
      for (replicate in seq_len(bootstrap_replicates)) {
        result <- bootstrap_qmm_once(
          fit,
          horizon = horizons[column],
          portmanteau_lag = portmanteau_lag,
          max_draw_attempts = max_draw_attempts
        )
        bootstrap[replicate, column] <- result[["statistic"]]
        attempts[replicate, column] <- result[["draw_attempts"]]
        selected_models[replicate, column] <- result[["selected_model"]]
      }
    }
  })

  limits <- vapply(
    seq_along(horizons),
    function(column) {
      as.numeric(stats::quantile(
        bootstrap[, column],
        probs = confidence_level,
        names = FALSE,
        type = 7L
      ))
    },
    numeric(1)
  )

  diagnostics <- data.frame(
    horizon = horizons,
    mean_draw_attempts = colMeans(attempts),
    rejected_draws = colSums(attempts - 1L),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      fit = fit,
      limits = data.frame(horizon = horizons, control_limit = limits),
      bootstrap = if (isTRUE(keep_bootstrap)) bootstrap else NULL,
      diagnostics = diagnostics,
      selected_models = if (isTRUE(keep_bootstrap)) selected_models else NULL,
      bootstrap_replicates = bootstrap_replicates,
      confidence_level = confidence_level,
      portmanteau_lag = portmanteau_lag,
      seed = seed,
      call = match.call()
    ),
    class = c("qmm_calibration", "fromage_chart")
  )
}
