#' Convert power-table output to plotting data
#'
#' @param resultados Output from `online_power_table_mm()`.
#' @param lags Numeric vector of lags used in the experiment.
#' @param h1_labels Labels for the phase 2 models.
#' @param tipo_teste Test type to extract.
#' @return A tidy data frame with power estimates and confidence intervals.
#' @export
cria_df_poder <- function(resultados,
                          lags,
                          h1_labels,
                          tipo_teste = c("Q Multi-Modelo", "Res\u00edduo Auto-Arima")) {
  tipo_teste <- match.arg(tipo_teste)

  if (tipo_teste == "Q Multi-Modelo") {
    medias <- resultados$resultados_q$media
    li <- resultados$resultados_q$li
    ls <- resultados$resultados_q$ls
  } else {
    medias <- resultados$resultados_res$media
    li <- resultados$resultados_res$li
    ls <- resultados$resultados_res$ls
  }

  if (length(h1_labels) != nrow(medias)) {
    stop("h1_labels must have one label for each row in resultados.", call. = FALSE)
  }

  df_media <- as.data.frame(medias)
  colnames(df_media) <- paste0("lag_", lags)
  df_media$modelo <- h1_labels

  df_li <- as.data.frame(li)
  colnames(df_li) <- paste0("lag_", lags)
  df_li$modelo <- h1_labels

  df_ls <- as.data.frame(ls)
  colnames(df_ls) <- paste0("lag_", lags)
  df_ls$modelo <- h1_labels

  df_media_long <- tidyr::pivot_longer(
    df_media,
    cols = dplyr::starts_with("lag_"),
    names_to = "Lag",
    values_to = "Poder"
  )

  df_li_long <- tidyr::pivot_longer(
    df_li,
    cols = dplyr::starts_with("lag_"),
    names_to = "Lag",
    values_to = "LI"
  )

  df_ls_long <- tidyr::pivot_longer(
    df_ls,
    cols = dplyr::starts_with("lag_"),
    names_to = "Lag",
    values_to = "LS"
  )

  df_final <- df_media_long
  df_final$Lag <- as.numeric(gsub("lag_", "", df_final$Lag))
  df_final$Tipo <- tipo_teste

  dplyr::bind_cols(
    df_final,
    LI = df_li_long$LI,
    LS = df_ls_long$LS
  )
}
