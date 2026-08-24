#' Validate a univariate numeric series
#'
#' Converts a numeric vector or univariate `ts` object to a plain numeric
#' vector after checking its length and confirming that every observation is
#' finite and non-missing. Matrices and multivariate time series are rejected.
#'
#' @param x Object to validate.
#' @param name Argument name used in error messages.
#' @param min_length Minimum permitted number of observations.
#'
#' @return A plain numeric vector.
#' @keywords internal
validate_numeric_series <- function(x, name, min_length = 2L) {
  is_multivariate <- !is.null(dim(x)) && !(inherits(x, "ts") && NCOL(x) == 1L)
  if (!is.numeric(x) || is_multivariate) {
    stop(name, " must be a numeric vector or a univariate time series.", call. = FALSE)
  }

  x <- as.numeric(x)
  if (length(x) < min_length) {
    stop(name, " must contain at least ", min_length, " observations.", call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop(name, " must contain only finite, non-missing values.", call. = FALSE)
  }

  x
}

#' Validate and convert one whole number
#'
#' Checks a scalar count or index before converting it to R's integer type.
#' Values above `.Machine$integer.max` are rejected to prevent conversion to
#' `NA`.
#'
#' @param x Value to validate.
#' @param name Argument name used in error messages.
#' @param minimum Smallest permitted value.
#'
#' @return A length-one integer.
#' @keywords internal
validate_whole_number <- function(x, name, minimum = 1L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x != floor(x)) {
    stop(name, " must be a single whole number.", call. = FALSE)
  }
  if (x < minimum) {
    stop(name, " must be at least ", minimum, ".", call. = FALSE)
  }
  if (x > .Machine$integer.max) {
    stop(name, " must not exceed ", .Machine$integer.max, ".", call. = FALSE)
  }

  as.integer(x)
}

#' Validate a probability
#'
#' Requires one finite numeric value strictly between zero and one.
#'
#' @param x Value to validate.
#' @param name Argument name used in error messages.
#'
#' @return The validated numeric value.
#' @keywords internal
validate_probability <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0 || x >= 1) {
    stop(name, " must be a single number strictly between 0 and 1.", call. = FALSE)
  }

  x
}

#' Validate a QMM Phase I fit object
#'
#' Checks the S3 class required by functions that consume [fit_qmm()] output.
#'
#' @param x Object to validate.
#'
#' @return `x` invisibly. An error is raised if `x` does not inherit from
#'   `qmm_fit`.
#' @keywords internal
validate_qmm_fit <- function(x) {
  if (!inherits(x, "qmm_fit")) {
    stop("fit must be an object created by fit_qmm().", call. = FALSE)
  }
  invisible(x)
}

#' Validate a QMM calibration object
#'
#' Checks the S3 class required by functions that consume [calibrate_qmm()]
#' output.
#'
#' @param x Object to validate.
#'
#' @return `x` invisibly. An error is raised if `x` does not inherit from
#'   `qmm_calibration`.
#' @keywords internal
validate_qmm_calibration <- function(x) {
  if (!inherits(x, "qmm_calibration")) {
    stop("calibration must be an object created by calibrate_qmm().", call. = FALSE)
  }
  invisible(x)
}

#' Evaluate code with a local reproducible seed
#'
#' When `seed` is supplied, temporarily sets R's random-number seed, evaluates
#' `code`, and restores the caller's previous random-number state even if the
#' expression fails. With `seed = NULL`, the expression uses and advances the
#' caller's random-number stream normally.
#'
#' @param seed Non-negative whole-number seed, or `NULL`.
#' @param code Expression to evaluate, supplied lazily.
#'
#' @return The value produced by `code`.
#' @keywords internal
with_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(force(code))
  }

  seed <- validate_whole_number(seed, "seed", minimum = 0L)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }

  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  force(code)
}

#' Evaluate an expression while collecting warnings
#'
#' Records unique warning messages and muffles them so a candidate-model fit
#' can preserve diagnostics without printing warnings during the fitting loop.
#' Errors are not caught by this helper.
#'
#' @param code Expression to evaluate, supplied lazily.
#'
#' @return A list with `value`, the expression result, and `warnings`, a
#'   character vector of unique warning messages.
#' @keywords internal
capture_warnings <- function(code) {
  warnings <- character()
  value <- withCallingHandlers(
    code,
    warning = function(warning) {
      warnings <<- c(warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )

  list(value = value, warnings = unique(warnings))
}

#' Check stationarity and invertibility of ARMA coefficients
#'
#' Uses the roots of the autoregressive and moving-average characteristic
#' polynomials. A coefficient vector is accepted only when every relevant root
#' lies strictly outside the unit circle by more than `tolerance`.
#'
#' @param coefficients Named ARMA coefficient vector.
#' @param p Autoregressive order.
#' @param q Moving-average order.
#' @param tolerance Positive separation required beyond modulus one.
#'
#' @return A length-one logical value. Missing required coefficient names yield
#'   `FALSE`.
#' @keywords internal
arma_is_admissible <- function(coefficients, p, q, tolerance = 1e-7) {
  ar <- if (p > 0L) coefficients[paste0("ar", seq_len(p))] else numeric()
  ma <- if (q > 0L) coefficients[paste0("ma", seq_len(q))] else numeric()

  if (anyNA(c(ar, ma))) {
    return(FALSE)
  }

  ar_ok <- p == 0L || all(Mod(polyroot(c(1, -ar))) > 1 + tolerance)
  ma_ok <- q == 0L || all(Mod(polyroot(c(1, ma))) > 1 + tolerance)

  ar_ok && ma_ok
}

#' Validate and order a named ARMA coefficient vector
#'
#' Requires finite coefficients named `ar1`, ..., `arp`, `ma1`, ..., `maq` as
#' appropriate for the requested orders. Extra named coefficients are dropped,
#' and the required values are returned in the canonical AR-then-MA order.
#'
#' @param coefficients Numeric coefficient vector.
#' @param p Autoregressive order.
#' @param q Moving-average order.
#' @param name Argument name used in error messages.
#'
#' @return A named numeric vector containing exactly the required coefficients.
#' @keywords internal
validate_coefficients <- function(coefficients, p, q, name = "coefficients") {
  required <- c(
    if (p > 0L) paste0("ar", seq_len(p)),
    if (q > 0L) paste0("ma", seq_len(q))
  )

  if (!is.numeric(coefficients) || anyNA(coefficients) || any(!is.finite(coefficients))) {
    stop(name, " must be a finite numeric vector.", call. = FALSE)
  }
  if (length(required) > 0L && is.null(names(coefficients))) {
    stop(name, " must use names such as ar1 and ma1.", call. = FALSE)
  }

  missing <- setdiff(required, names(coefficients))
  if (length(missing) > 0L) {
    stop(name, " is missing: ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }

  coefficients[required]
}
