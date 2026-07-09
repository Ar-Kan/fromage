#' Build a grid of ARMA candidate models
#'
#' @param max.p Maximum AR order.
#' @param max.q Maximum MA order.
#' @param include_white_noise Whether to include ARMA(0, 0).
#' @return A data frame with columns `p` and `q`.
#' @export
arma_candidates <- function(max.p = 2, max.q = 2, include_white_noise = FALSE) {
  validate_count(max.p, "max.p", allow_zero = TRUE)
  validate_count(max.q, "max.q", allow_zero = TRUE)

  modelos_candidatos <- expand.grid(
    p = 0:max.p,
    q = 0:max.q
  )

  if (!include_white_noise) {
    modelos_candidatos <- modelos_candidatos[
      !(modelos_candidatos$p == 0 & modelos_candidatos$q == 0),
      ,
      drop = FALSE
    ]
  }

  rownames(modelos_candidatos) <- NULL
  modelos_candidatos
}

#' Fit a multimodel ARMA candidate set
#'
#' @param serie Numeric time series.
#' @param modelos_df Data frame with candidate columns `p` and `q`.
#' @param criterio Information criterion stored in fitted `forecast::Arima`
#'   objects, usually `"aic"`, `"aicc"`, or `"bic"`.
#' @param cutoff Optional delta-IC cutoff for retaining substantial models.
#' @return A list containing fitted models, model weights, IC values, and the
#'   candidate grid used.
#' @export
ajustar_modelos <- function(serie,
                            modelos_df,
                            criterio = "aic",
                            cutoff = NULL) {
  if (!is.numeric(serie) || anyNA(serie)) {
    stop("serie must be a numeric vector without missing values.", call. = FALSE)
  }
  if (!is.data.frame(modelos_df) || !all(c("p", "q") %in% names(modelos_df))) {
    stop("modelos_df must be a data frame with columns p and q.", call. = FALSE)
  }
  if (!is.character(criterio) || length(criterio) != 1) {
    stop("criterio must be a single character value.", call. = FALSE)
  }
  if (!is.null(cutoff) && (!is.numeric(cutoff) || length(cutoff) != 1 || cutoff < 0)) {
    stop("cutoff must be NULL or a non-negative number.", call. = FALSE)
  }

  resultados <- vector("list", nrow(modelos_df))
  ics <- rep(Inf, nrow(modelos_df))

  for (i in seq_len(nrow(modelos_df))) {
    p <- modelos_df$p[i]
    q <- modelos_df$q[i]

    tryCatch({
      modelo <- forecast::Arima(
        serie,
        order = c(p, 0, q),
        include.mean = FALSE
      )

      if (is.null(modelo[[criterio]])) {
        stop("criterion not found")
      }

      resultados[[i]] <- modelo
      ics[i] <- modelo[[criterio]]
    }, error = function(e) {
      resultados[[i]] <<- NULL
      ics[i] <<- Inf
    })
  }

  if (!any(is.finite(ics))) {
    stop("No candidate model could be fitted.", call. = FALSE)
  }

  delta_ic <- ics - min(ics, na.rm = TRUE)
  pesos <- exp(-delta_ic / 2) / sum(exp(-delta_ic / 2), na.rm = TRUE)

  if (is.null(cutoff)) {
    return(list(
      modelos = resultados,
      pesos = pesos,
      ics = ics,
      modelos_df = modelos_df
    ))
  }

  substantial <- which(delta_ic <= cutoff)
  pesos_cut <- exp(-(ics[substantial] - min(ics[substantial])) / 2)
  pesos_cut <- pesos_cut / sum(pesos_cut)

  list(
    modelos = resultados[substantial],
    pesos = pesos_cut,
    ics = ics[substantial],
    modelos_df = modelos_df[substantial, , drop = FALSE]
  )
}
