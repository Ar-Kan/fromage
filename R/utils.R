#' Validate a single positive whole number
#'
#' @param x Value to validate.
#' @param name Argument name used in error messages.
#' @param allow_zero Whether zero is accepted.
#' @return The validated value, invisibly.
#' @noRd
validate_count <- function(x, name, allow_zero = FALSE) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x != floor(x)) {
    stop(name, " must be a single whole number.", call. = FALSE)
  }

  if (allow_zero) {
    valid <- x >= 0
  } else {
    valid <- x > 0
  }

  if (!valid) {
    stop(name, " is outside the allowed range.", call. = FALSE)
  }

  invisible(x)
}

#' Validate an ARMA order vector
#'
#' @param order Numeric vector of length three: `c(p, d, q)`.
#' @return The validated order.
#' @noRd
validate_arma_order <- function(order) {
  if (!is.numeric(order) || length(order) != 3 || anyNA(order)) {
    stop("order must be a numeric vector of length 3: c(p, d, q).", call. = FALSE)
  }

  if (any(order != floor(order)) || any(order < 0)) {
    stop("order values must be non-negative whole numbers.", call. = FALSE)
  }

  order
}

#' Validate ARMA coefficient names
#'
#' @param coefs Named numeric vector.
#' @param p AR order.
#' @param q MA order.
#' @param name Argument name used in error messages.
#' @return The validated coefficient vector.
#' @noRd
validate_arma_coefs <- function(coefs, p, q, name = "coefs") {
  if (!is.numeric(coefs) || is.null(names(coefs)) || anyNA(coefs)) {
    stop(name, " must be a named numeric vector without missing values.", call. = FALSE)
  }

  required <- c(
    if (p > 0) paste0("ar", seq_len(p)),
    if (q > 0) paste0("ma", seq_len(q))
  )

  missing <- setdiff(required, names(coefs))
  if (length(missing) > 0) {
    stop(
      name,
      " is missing required coefficient(s): ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  coefs
}

#' Calculate Wilson confidence intervals for binary violations
#'
#' @param violations Numeric or logical vector of binary violations.
#' @return A list with lower bound, mean, and upper bound.
#' @noRd
wilson_summary <- function(violations) {
  violations <- as.numeric(violations)
  ic <- binom::binom.confint(
    sum(violations),
    length(violations),
    methods = "wilson"
  )

  list(
    li = ic$lower,
    media = ic$mean,
    ls = ic$upper
  )
}
