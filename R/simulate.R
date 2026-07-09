#' Simulate future observations from an ARMA model conditional on history
#'
#' @param n Number of observations to simulate.
#' @param coefs Named ARMA coefficient vector. Names must include `ar1`, ...,
#'   `arp` and/or `ma1`, ..., `maq` according to `ordem`.
#' @param ordem ARIMA order vector `c(p, d, q)`. Only `p` and `q` are used.
#' @param y Historical observations used to condition the simulation.
#' @param eps Historical innovations used to condition the simulation.
#' @param sd Innovation standard deviation.
#' @return Numeric vector with `n` simulated observations.
#' @export
arima_structbreak_sim <- function(n,
                                  coefs,
                                  ordem,
                                  y,
                                  eps,
                                  sd = 1) {
  validate_count(n, "n")
  ordem <- validate_arma_order(ordem)
  p <- ordem[1]
  q <- ordem[3]
  validate_arma_coefs(coefs, p, q)

  if (!is.numeric(y) || length(y) < p || anyNA(y)) {
    stop("y must be a numeric history with enough observations for the AR order.", call. = FALSE)
  }
  if (!is.numeric(eps) || length(eps) < q || anyNA(eps)) {
    stop("eps must be a numeric history with enough observations for the MA order.", call. = FALSE)
  }
  if (!is.numeric(sd) || length(sd) != 1 || is.na(sd) || sd <= 0) {
    stop("sd must be a positive number.", call. = FALSE)
  }

  repeat {
    result <- tryCatch({
      ar_coefs <- if (p > 0) coefs[paste0("ar", seq_len(p))] else numeric(0)
      ma_coefs <- if (q > 0) coefs[paste0("ma", seq_len(q))] else numeric(0)

      novos_eps <- stats::rnorm(n, sd = sd)
      novo_y <- numeric(n)

      y_hist <- as.numeric(y)
      eps_hist <- as.numeric(eps)

      for (t in seq_len(n)) {
        ar_part <- 0
        if (p > 0) {
          for (i in seq_len(p)) {
            ar_part <- ar_part + ar_coefs[i] * y_hist[length(y_hist) - i + 1]
          }
        }

        ma_part <- 0
        if (q > 0) {
          for (j in seq_len(q)) {
            ma_part <- ma_part + ma_coefs[j] * eps_hist[length(eps_hist) - j + 1]
          }
        }

        novo_y[t] <- ar_part + novos_eps[t] + ma_part
        y_hist <- c(y_hist, novo_y[t])
        eps_hist <- c(eps_hist, novos_eps[t])
      }

      novo_y
    }, error = function(e) NULL)

    if (!is.null(result)) {
      break
    }
  }

  result
}

#' Simulate an ARMA process with one structural break
#'
#' @param n Number of observations in each phase.
#' @param fase1_coefs Named coefficient vector for phase 1.
#' @param fase2_coefs Named coefficient vector for phase 2.
#' @param order ARIMA order vector `c(p, d, q)`. Only `p` and `q` are used.
#' @param sd Innovation standard deviation.
#' @param burn_in Number of observations discarded before phase 1.
#' @return A list with phase series, phase innovations, and an `auto.arima` fit
#'   estimated on phase 1.
#' @export
simula_processo <- function(n,
                            fase1_coefs,
                            fase2_coefs,
                            order = c(1, 0, 1),
                            sd = 1,
                            burn_in = 1000) {
  validate_count(n, "n")
  validate_count(burn_in, "burn_in", allow_zero = TRUE)
  order <- validate_arma_order(order)
  p <- order[1]
  q <- order[3]
  validate_arma_coefs(fase1_coefs, p, q, "fase1_coefs")
  validate_arma_coefs(fase2_coefs, p, q, "fase2_coefs")

  if (!is.numeric(sd) || length(sd) != 1 || is.na(sd) || sd <= 0) {
    stop("sd must be a positive number.", call. = FALSE)
  }

  repeat {
    result <- tryCatch({
      n_total <- burn_in + 2 * n
      eps <- stats::rnorm(n_total, sd = sd)
      y <- numeric(n_total)

      start_t <- max(p, q) + 1
      for (t in start_t:n_total) {
        if (t <= burn_in + n) {
          ar_coefs <- if (p > 0) fase1_coefs[paste0("ar", seq_len(p))] else numeric(0)
          ma_coefs <- if (q > 0) fase1_coefs[paste0("ma", seq_len(q))] else numeric(0)
        } else {
          ar_coefs <- if (p > 0) fase2_coefs[paste0("ar", seq_len(p))] else numeric(0)
          ma_coefs <- if (q > 0) fase2_coefs[paste0("ma", seq_len(q))] else numeric(0)
        }

        ar_part <- 0
        if (p > 0) {
          for (i in seq_len(p)) {
            ar_part <- ar_part + ar_coefs[i] * y[t - i]
          }
        }

        ma_part <- 0
        if (q > 0) {
          for (j in seq_len(q)) {
            ma_part <- ma_part + ma_coefs[j] * eps[t - j]
          }
        }

        y[t] <- ar_part + eps[t] + ma_part
      }

      y_fase1 <- y[(burn_in + 1):(burn_in + n)]
      y_fase2 <- y[(burn_in + n + 1):(burn_in + 2 * n)]
      eps_fase1 <- eps[(burn_in + 1):(burn_in + n)]
      eps_fase2 <- eps[(burn_in + n + 1):(burn_in + 2 * n)]

      fit <- forecast::auto.arima(
        y_fase1,
        seasonal = FALSE,
        stationary = TRUE,
        allowmean = FALSE,
        allowdrift = FALSE
      )

      list(
        fase1 = y_fase1,
        fase2 = y_fase2,
        eps_fase1 = eps_fase1,
        eps_fase2 = eps_fase2,
        fit = fit
      )
    }, error = function(e) NULL)

    if (!is.null(result)) {
      break
    }
  }

  result
}
