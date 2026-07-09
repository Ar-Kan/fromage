#' Calculate a multimodel Ljung-Box Q statistic
#'
#' @param serie Numeric series to evaluate.
#' @param resultados_mmi Output from `ajustar_modelos()`.
#' @param lag_q Lag used in `stats::Box.test()`.
#' @return Numeric multimodel Q statistic.
#' @export
multimodel_Q <- function(serie,
                         resultados_mmi,
                         lag_q = 10) {
  if (!is.numeric(serie) || anyNA(serie)) {
    stop("serie must be a numeric vector without missing values.", call. = FALSE)
  }
  validate_count(lag_q, "lag_q")

  modelos_validos <- which(resultados_mmi$pesos > 0)
  pesos <- resultados_mmi$pesos[modelos_validos]
  stats <- numeric(length(modelos_validos))

  for (i in seq_along(modelos_validos)) {
    m <- modelos_validos[i]
    modelo <- resultados_mmi$modelos[[m]]

    ajuste <- forecast::Arima(
      serie,
      include.mean = FALSE,
      order = forecast::arimaorder(modelo),
      fixed = modelo$coef
    )

    q_i <- stats::Box.test(
      ajuste$residuals,
      lag = lag_q,
      type = "Ljung-Box"
    )$statistic

    stats[i] <- q_i
  }

  pesos <- pesos / sum(pesos)
  as.numeric(sum(stats * pesos))
}

#' Bootstrap the multimodel Q statistic
#'
#' @param resultados_mmi Output from `ajustar_modelos()`.
#' @param fase1 Phase 1 observations.
#' @param eps_fase1 Phase 1 innovations.
#' @param lag Number of observations replaced by simulated values.
#' @param n_boot Number of bootstrap samples.
#' @param lag_q Lag used in the Ljung-Box statistic.
#' @return Numeric vector of bootstrap Q statistics.
#' @export
bootstrap_multimodelo <- function(resultados_mmi,
                                  fase1,
                                  eps_fase1,
                                  lag,
                                  n_boot = 500,
                                  lag_q = 10) {
  validate_count(lag, "lag")
  validate_count(n_boot, "n_boot")
  validate_count(lag_q, "lag_q")

  if (!is.numeric(fase1) || anyNA(fase1)) {
    stop("fase1 must be a numeric vector without missing values.", call. = FALSE)
  }
  if (!is.numeric(eps_fase1) || anyNA(eps_fase1)) {
    stop("eps_fase1 must be a numeric vector without missing values.", call. = FALSE)
  }

  pesos <- resultados_mmi$pesos
  modelos_validos <- which(pesos > 0)
  n <- length(fase1)
  boot <- numeric(n_boot)

  for (b in seq_len(n_boot)) {
    modelo_idx <- sample(
      modelos_validos,
      size = 1,
      prob = pesos[modelos_validos]
    )

    modelo <- resultados_mmi$modelos[[modelo_idx]]
    ordem <- forecast::arimaorder(modelo)
    success <- FALSE
    attempt <- 0

    while (!success && attempt < 1000) {
      attempt <- attempt + 1

      tryCatch({
        if (length(modelo$coef) == 1) {
          coef_star <- stats::rnorm(
            1,
            modelo$coef,
            sqrt(modelo$var.coef)
          )
          names(coef_star) <- names(modelo$coef)
        } else {
          coef_star <- mvtnorm::rmvnorm(
            1,
            modelo$coef,
            modelo$var.coef
          )[1, ]
        }

        y_star <- arima_structbreak_sim(
          n = lag,
          coefs = coef_star,
          ordem = ordem,
          y = fase1,
          eps = eps_fase1,
          sd = sqrt(modelo$sigma2)
        )

        success <- TRUE
      }, error = function(e) {})
    }

    if (!success) {
      stop("Bootstrap failed after 1000 attempts.", call. = FALSE)
    }

    if (lag < n) {
      colagem <- c(fase1[(lag + 1):n], y_star)
    } else {
      colagem <- y_star
    }

    boot[b] <- multimodel_Q(
      serie = colagem,
      resultados_mmi = resultados_mmi,
      lag_q = lag_q
    )
  }

  boot
}
