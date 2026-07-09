#' Run a Monte Carlo power test for multimodel Q and residual tests
#'
#' @inheritParams parallel_power_test
#' @return A list with binary violation vectors for each test.
#' @export
power_test <- function(n,
                       lag,
                       n_mc,
                       n_boot = 500,
                       alpha = 0.05,
                       lag_q = 10,
                       max.p = 2,
                       max.q = 2,
                       criterio = "aicc",
                       cutoff = NULL,
                       order = c(1, 0, 1),
                       h0_model,
                       h1_model,
                       burn_in = 1000) {
  validate_count(n, "n")
  validate_count(lag, "lag")
  validate_count(n_mc, "n_mc")
  validate_count(n_boot, "n_boot")
  validate_count(lag_q, "lag_q")
  validate_count(max.p, "max.p", allow_zero = TRUE)
  validate_count(max.q, "max.q", allow_zero = TRUE)

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("alpha must be between 0 and 1.", call. = FALSE)
  }

  order <- validate_arma_order(order)
  validate_arma_coefs(h0_model, order[1], order[3], "h0_model")
  validate_arma_coefs(h1_model, order[1], order[3], "h1_model")

  violations_q <- numeric(n_mc)
  violations_res <- numeric(n_mc)
  modelos_candidatos <- arma_candidates(max.p, max.q)

  for (m in seq_len(n_mc)) {
    processo <- simula_processo(
      n = n,
      fase1_coefs = h0_model,
      fase2_coefs = h1_model,
      order = order,
      burn_in = burn_in
    )

    h0 <- processo$fase1
    h0_eps <- processo$eps_fase1

    h0_mmi <- ajustar_modelos(
      serie = h0,
      modelos_df = modelos_candidatos,
      criterio = criterio,
      cutoff = cutoff
    )

    fit_auto <- forecast::auto.arima(
      h0,
      seasonal = FALSE,
      stationary = TRUE,
      allowmean = FALSE,
      allowdrift = FALSE
    )

    if (lag < n) {
      h1 <- processo$fase2[seq_len(lag)]
      colagem <- c(h0[(lag + 1):n], h1)
    } else {
      h1 <- processo$fase2
      colagem <- h1
    }

    q_obs <- multimodel_Q(
      serie = colagem,
      resultados_mmi = h0_mmi,
      lag_q = lag_q
    )

    boot_q <- bootstrap_multimodelo(
      resultados_mmi = h0_mmi,
      fase1 = h0,
      eps_fase1 = h0_eps,
      lag = lag,
      n_boot = n_boot,
      lag_q = lag_q
    )

    ul_q <- stats::quantile(
      boot_q,
      1 - alpha,
      na.rm = TRUE
    )

    violations_q[m] <- q_obs > ul_q

    ajuste_res <- forecast::Arima(
      colagem,
      include.mean = FALSE,
      order = forecast::arimaorder(fit_auto),
      fixed = fit_auto$coef
    )

    residuo_final <- as.numeric(utils::tail(ajuste_res$residuals, 1))

    violations_res[m] <-
      residuo_final < stats::qnorm(alpha / 2) |
        residuo_final > stats::qnorm(1 - alpha / 2)
  }

  list(
    violacoes_q = violations_q,
    violacoes_res = violations_res
  )
}

#' Run a parallel Monte Carlo power test
#'
#' @param n Number of observations in each phase.
#' @param lag Number of phase 2 observations introduced into the monitored
#'   series.
#' @param n_mc Number of Monte Carlo replications.
#' @param n_boot Number of bootstrap replications per Monte Carlo run.
#' @param alpha Test significance level.
#' @param lag_q Lag used in the Ljung-Box statistic.
#' @param max.p Maximum AR order in the candidate set.
#' @param max.q Maximum MA order in the candidate set.
#' @param criterio Information criterion used to weight models.
#' @param cutoff Optional delta-IC cutoff for substantial models.
#' @param order True ARIMA order vector used in simulation.
#' @param h0_model Named coefficient vector for phase 1.
#' @param h1_model Named coefficient vector for phase 2.
#' @param burn_in Number of burn-in observations.
#' @param cores Number of worker processes. Use `1` for sequential execution.
#' @param progress Whether to show a text progress bar in parallel execution.
#' @return A list with Wilson confidence summaries for both tests.
#' @importFrom foreach %dopar%
#' @export
parallel_power_test <- function(n,
                                lag,
                                n_mc,
                                n_boot = 500,
                                alpha = 0.05,
                                lag_q = 10,
                                max.p = 2,
                                max.q = 2,
                                criterio = "aicc",
                                cutoff = NULL,
                                order = c(1, 0, 1),
                                h0_model,
                                h1_model,
                                burn_in = 1000,
                                cores = max(1, parallel::detectCores() - 1),
                                progress = interactive()) {
  validate_count(cores, "cores")

  if (cores == 1) {
    teste <- power_test(
      n = n,
      lag = lag,
      n_mc = n_mc,
      n_boot = n_boot,
      alpha = alpha,
      lag_q = lag_q,
      max.p = max.p,
      max.q = max.q,
      criterio = criterio,
      cutoff = cutoff,
      order = order,
      h0_model = h0_model,
      h1_model = h1_model,
      burn_in = burn_in
    )
  } else {
    cl <- parallel::makeCluster(cores, type = "SOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)

    pb <- NULL
    if (isTRUE(progress)) {
      pb <- utils::txtProgressBar(max = n_mc, style = 3)
      on.exit(close(pb), add = TRUE)
    }

    progress_fun <- function(n) {
      if (!is.null(pb)) {
        utils::setTxtProgressBar(pb, n)
      }
    }

    doSNOW::registerDoSNOW(cl)
    opts <- if (isTRUE(progress)) list(progress = progress_fun) else list()

    teste_raw <- foreach::foreach(
      i = seq_len(n_mc),
      .combine = rbind,
      .options.snow = opts,
      .export = c(
        "validate_count",
        "validate_arma_order",
        "validate_arma_coefs",
        "wilson_summary",
        "arima_structbreak_sim",
        "simula_processo",
        "arma_candidates",
        "ajustar_modelos",
        "multimodel_Q",
        "bootstrap_multimodelo",
        "power_test"
      ),
      .packages = c(
        "binom",
        "forecast",
        "mvtnorm"
      )
    ) %dopar% {
      power_test(
        n = n,
        lag = lag,
        n_mc = 1,
        n_boot = n_boot,
        alpha = alpha,
        lag_q = lag_q,
        max.p = max.p,
        max.q = max.q,
        criterio = criterio,
        cutoff = cutoff,
        order = order,
        h0_model = h0_model,
        h1_model = h1_model,
        burn_in = burn_in
      )
    }

    teste <- list(
      violacoes_q = as.numeric(teste_raw[, "violacoes_q"]),
      violacoes_res = as.numeric(teste_raw[, "violacoes_res"])
    )
  }

  list(
    resultados_q = wilson_summary(teste$violacoes_q),
    resultados_res = wilson_summary(teste$violacoes_res)
  )
}

#' Build an online power table across lags and phase 2 models
#'
#' @inheritParams parallel_power_test
#' @param lags Numeric vector of lag values.
#' @param h1_models List of named phase 2 coefficient vectors.
#' @return A nested list with matrices of Wilson lower bounds, means, and upper
#'   bounds for both tests.
#' @export
online_power_table_mm <- function(lags,
                                  h0_model,
                                  h1_models,
                                  n,
                                  n_mc,
                                  n_boot = 500,
                                  alpha = 0.05,
                                  lag_q = 10,
                                  max.p = 2,
                                  max.q = 2,
                                  criterio = "aicc",
                                  cutoff = NULL,
                                  order = c(1, 0, 1),
                                  burn_in = 1000,
                                  cores = max(1, parallel::detectCores() - 1),
                                  progress = interactive()) {
  if (!is.numeric(lags) || anyNA(lags) || any(lags <= 0)) {
    stop("lags must be a positive numeric vector.", call. = FALSE)
  }
  if (!is.list(h1_models) || length(h1_models) == 0) {
    stop("h1_models must be a non-empty list.", call. = FALSE)
  }

  matriz_medias_q <- matrix(0, nrow = length(h1_models), ncol = length(lags))
  matriz_li_q <- matrix(0, nrow = length(h1_models), ncol = length(lags))
  matriz_ls_q <- matrix(0, nrow = length(h1_models), ncol = length(lags))

  matriz_medias_res <- matrix(0, nrow = length(h1_models), ncol = length(lags))
  matriz_li_res <- matrix(0, nrow = length(h1_models), ncol = length(lags))
  matriz_ls_res <- matrix(0, nrow = length(h1_models), ncol = length(lags))

  for (i in seq_along(h1_models)) {
    for (j in seq_along(lags)) {
      resultado_ij <- parallel_power_test(
        n = n,
        lag = lags[j],
        n_mc = n_mc,
        n_boot = n_boot,
        alpha = alpha,
        lag_q = lag_q,
        max.p = max.p,
        max.q = max.q,
        criterio = criterio,
        cutoff = cutoff,
        order = order,
        h0_model = h0_model,
        h1_model = h1_models[[i]],
        burn_in = burn_in,
        cores = cores,
        progress = progress
      )

      matriz_medias_q[i, j] <- resultado_ij$resultados_q$media
      matriz_li_q[i, j] <- resultado_ij$resultados_q$li
      matriz_ls_q[i, j] <- resultado_ij$resultados_q$ls

      matriz_medias_res[i, j] <- resultado_ij$resultados_res$media
      matriz_li_res[i, j] <- resultado_ij$resultados_res$li
      matriz_ls_res[i, j] <- resultado_ij$resultados_res$ls
    }
  }

  list(
    resultados_q = list(
      li = matriz_li_q,
      media = matriz_medias_q,
      ls = matriz_ls_q
    ),
    resultados_res = list(
      li = matriz_li_res,
      media = matriz_medias_res,
      ls = matriz_ls_res
    )
  )
}
