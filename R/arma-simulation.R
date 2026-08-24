#' Simulate an ARMA process with a coefficient change
#'
#' Simulates one continuous series. The first segment uses the reference ARMA
#' coefficients, and the second segment continues from the same history using
#' the monitoring coefficients.
#'
#' @param reference_length Number of returned reference observations (Phase I).
#' @param monitoring_length Number of returned monitoring observations (Phase II).
#' @param reference_coefficients Named ARMA coefficients for the reference
#'   segment, such as `c(ar1 = 0.4, ma1 = 0.3)`.
#' @param monitoring_coefficients Named ARMA coefficients for the monitoring
#'   segment.
#' @param order Integer vector `c(p, q)`.
#' @param innovation_sd Innovation standard deviation.
#' @param burn_in Number of initial observations discarded.
#' @param seed Optional random-number seed.
#'
#' @return A list containing `series`, the complete continuous simulation;
#'   `reference` and `monitoring`, its Phase I and Phase II segments;
#'   `reference_innovations` and `monitoring_innovations`, the corresponding
#'   simulated innovations; and `order`, the named ARMA order. The innovations
#'   are returned for simulation validation and are not required by the public
#'   QMM workflow.
#' @export
simulate_arma_change <- function(reference_length,
                                 monitoring_length,
                                 reference_coefficients,
                                 monitoring_coefficients,
                                 order = c(1L, 1L),
                                 innovation_sd = 1,
                                 burn_in = 1000L,
                                 seed = NULL) {
  reference_length <- validate_whole_number(reference_length, "reference_length")
  monitoring_length <- validate_whole_number(monitoring_length, "monitoring_length")
  burn_in <- validate_whole_number(burn_in, "burn_in", minimum = 0L)
  if (!is.numeric(order) || length(order) != 2L || anyNA(order) ||
      any(order < 0) || any(order != floor(order))) {
    stop("order must be c(p, q) with non-negative whole numbers.", call. = FALSE)
  }
  p <- as.integer(order[1])
  q <- as.integer(order[2])
  reference_coefficients <- validate_coefficients(
    reference_coefficients,
    p,
    q,
    "reference_coefficients"
  )
  monitoring_coefficients <- validate_coefficients(
    monitoring_coefficients,
    p,
    q,
    "monitoring_coefficients"
  )
  if (!arma_is_admissible(reference_coefficients, p, q)) {
    stop("reference_coefficients must define a stationary and invertible model.", call. = FALSE)
  }
  if (!arma_is_admissible(monitoring_coefficients, p, q)) {
    stop("monitoring_coefficients must define a stationary and invertible model.", call. = FALSE)
  }
  if (!is.numeric(innovation_sd) || length(innovation_sd) != 1L ||
      is.na(innovation_sd) || !is.finite(innovation_sd) || innovation_sd <= 0) {
    stop("innovation_sd must be a positive finite number.", call. = FALSE)
  }

  total <- burn_in + reference_length + monitoring_length
  result <- with_seed(seed, {
    innovations <- stats::rnorm(total, sd = innovation_sd)
    observations <- numeric(total)

    for (time in seq_len(total)) {
      coefficients <- if (time <= burn_in + reference_length) {
        reference_coefficients
      } else {
        monitoring_coefficients
      }
      ar <- if (p > 0L) coefficients[paste0("ar", seq_len(p))] else numeric()
      ma <- if (q > 0L) coefficients[paste0("ma", seq_len(q))] else numeric()
      available_ar <- if (p > 0L) seq_len(min(p, time - 1L)) else integer()
      available_ma <- if (q > 0L) seq_len(min(q, time - 1L)) else integer()
      ar_part <- if (length(available_ar) > 0L) {
        sum(ar[available_ar] * observations[time - available_ar])
      } else {
        0
      }
      ma_part <- if (length(available_ma) > 0L) {
        sum(ma[available_ma] * innovations[time - available_ma])
      } else {
        0
      }
      observations[time] <- ar_part + ma_part + innovations[time]
    }

    keep <- (burn_in + 1L):total
    list(observations = observations[keep], innovations = innovations[keep])
  })

  reference_index <- seq_len(reference_length)
  monitoring_index <- reference_length + seq_len(monitoring_length)
  list(
    series = result$observations,
    reference = result$observations[reference_index],
    monitoring = result$observations[monitoring_index],
    reference_innovations = result$innovations[reference_index],
    monitoring_innovations = result$innovations[monitoring_index],
    order = c(p = p, q = q)
  )
}
