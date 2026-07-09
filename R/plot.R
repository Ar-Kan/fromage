#' Plot power curves for multimodel Q experiments
#'
#' @param df_plot Data frame produced by `cria_df_poder()` or a bound set of
#'   outputs from that function.
#' @param info_experimento Optional caption describing the experiment.
#' @param h1_labels Optional model order for the legend.
#' @param alpha_line Reference significance line.
#' @return A `ggplot` object.
#' @export
plot_power_curve <- function(df_plot,
                             info_experimento = NULL,
                             h1_labels = NULL,
                             alpha_line = 0.05) {
  required <- c("Lag", "Poder", "modelo", "LI", "LS", "Tipo")
  if (!is.data.frame(df_plot) || !all(required %in% names(df_plot))) {
    stop("df_plot must contain columns: ", paste(required, collapse = ", "), ".", call. = FALSE)
  }

  if (!is.null(h1_labels)) {
    df_plot$modelo <- factor(df_plot$modelo, levels = h1_labels)
  }

  ggplot2::ggplot(
    df_plot,
    ggplot2::aes(
      x = Lag,
      y = Poder,
      color = modelo,
      group = modelo
    )
  ) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = LI, ymax = LS),
      width = 2,
      alpha = 0.45
    ) +
    ggplot2::geom_hline(
      yintercept = alpha_line,
      linetype = "dashed",
      linewidth = 0.8
    ) +
    ggplot2::facet_wrap(~Tipo) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.1)
    ) +
    ggplot2::labs(
      x = "Lag",
      y = "Poder",
      color = "Modelo Fase 2",
      caption = info_experimento
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "gray90"),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      plot.caption = ggplot2::element_text(
        hjust = 0.5,
        size = 7,
        margin = ggplot2::margin(t = 15)
      )
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        nrow = 2,
        byrow = TRUE
      )
    )
}
