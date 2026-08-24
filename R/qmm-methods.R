#' Print a QMM Phase I fit
#'
#' Displays the information criterion, candidate counts, and the differences
#' and weights of retained models.
#'
#' @param x A `qmm_fit` object.
#' @param ... Additional arguments passed to [print()].
#'
#' @return `x` invisibly.
#' @method print qmm_fit
#' @export
print.qmm_fit <- function(x, ...) {
  cat("Multimodel Q reference fit\n")
  cat("Reference observations:", length(x$reference), "\n")
  cat("Criterion:", x$criterion, "\n")
  cat(
    "Candidates:", nrow(x$candidates),
    sprintf("(%d included)\n", sum(x$candidates$included))
  )
  visible <- x$candidates[x$candidates$included, c("model", "delta", "weight")]
  print(visible, row.names = FALSE, digits = 4)
  invisible(x)
}

#' Summarize a QMM Phase I fit
#'
#' @param object A `qmm_fit` object.
#' @param ... Unused.
#'
#' @return The complete candidate-model data frame, including failed and
#'   excluded candidates.
#' @method summary qmm_fit
#' @export
summary.qmm_fit <- function(object, ...) {
  object$candidates
}

#' Print a QMM bootstrap calibration
#'
#' Displays the number of bootstrap replicates, pointwise control-limit
#' probability, and the horizon-specific limits.
#'
#' @param x A `qmm_calibration` object.
#' @param ... Additional arguments passed to [print()].
#'
#' @return `x` invisibly.
#' @method print qmm_calibration
#' @export
print.qmm_calibration <- function(x, ...) {
  cat("Multimodel Q bootstrap calibration\n")
  cat("Bootstrap replicates per horizon:", x$bootstrap_replicates, "\n")
  cat("Pointwise control-limit probability:", x$confidence_level, "\n")
  print(x$limits, row.names = FALSE, digits = 5)
  invisible(x)
}

#' Summarize a QMM bootstrap calibration
#'
#' @param object A `qmm_calibration` object.
#' @param ... Unused.
#'
#' @return A data frame joining each horizon's control limit with its
#'   coefficient-draw diagnostics.
#' @method summary qmm_calibration
#' @export
summary.qmm_calibration <- function(object, ...) {
  merge(object$limits, object$diagnostics, by = "horizon", sort = FALSE)
}

#' Print a QMM monitoring result
#'
#' Displays the Phase I and Phase II sample sizes, pointwise control-limit
#' probability, signal count, and first signal horizon when a signal exists.
#'
#' @param x A `qmm_monitor` object.
#' @param ... Unused.
#'
#' @return `x` invisibly.
#' @method print qmm_monitor
#' @export
print.qmm_monitor <- function(x, ...) {
  signal_count <- sum(x$results$signal)
  cat("Multimodel Q (QMM) control chart\n")
  cat("Reference observations (Phase I):", x$reference_length, "\n")
  cat("Monitoring observations (Phase II):", nrow(x$results), "\n")
  cat(
    "Pointwise control-limit probability:",
    x$calibration$confidence_level,
    "\n"
  )
  cat("Phase II signals:", signal_count, "\n")
  if (signal_count > 0L) {
    cat("First signal at horizon:", min(x$results$horizon[x$results$signal]), "\n")
  }
  invisible(x)
}

#' Summarize a QMM monitoring result
#'
#' @param object A `qmm_monitor` object.
#' @param ... Unused.
#'
#' @return A data frame with the Phase II horizon, observation, QMM statistic,
#'   pointwise upper control limit, and signal indicator.
#' @method summary qmm_monitor
#' @export
summary.qmm_monitor <- function(object, ...) {
  object$results
}

#' Plot a multimodel Q monitoring result
#'
#' @param x A monitoring result created by [monitor_qmm()] or [qmm_chart()].
#' @param type Plot the Phase II QMM statistic and control limit
#'   (`"statistic"`), the complete Phase I/Phase II series (`"series"`), or an
#'   aligned two-panel view containing both (`"overview"`).
#' @param ... Unused.
#'
#' @return A `ggplot` object.
#' @export
plot.qmm_monitor <- function(x, type = c("statistic", "series", "overview"), ...) {
  type <- match.arg(type)

  if (type == "statistic") {
    plot_data <- x$results
    return(
      ggplot2::ggplot(plot_data, ggplot2::aes(x = horizon)) +
        ggplot2::geom_line(ggplot2::aes(y = control_limit), linetype = "dashed") +
        ggplot2::geom_line(ggplot2::aes(y = statistic), color = "#225EA8") +
        ggplot2::geom_point(
          ggplot2::aes(y = statistic, color = signal),
          size = 2.2
        ) +
        ggplot2::scale_color_manual(values = c(`FALSE` = "#225EA8", `TRUE` = "#D7301F")) +
        ggplot2::labs(
          x = "Monitoring horizon",
          y = "QMM statistic",
          color = "Signal",
          title = "Multimodel Q (QMM) control chart",
          subtitle = "Phase II; the dashed line is the pointwise upper control limit"
        ) +
        ggplot2::theme_minimal(base_size = 12)
    )
  }

  reference <- x$calibration$fit$reference
  series <- c(reference, x$new_data)
  reference_length <- length(reference)
  signal_positions <- reference_length + x$results$horizon[x$results$signal]
  plot_data <- data.frame(
    time = seq_along(series),
    value = series,
    signal = seq_along(series) %in% signal_positions
  )

  if (type == "series") {
    return(
      ggplot2::ggplot(plot_data, ggplot2::aes(x = time, y = value)) +
        ggplot2::geom_line(color = "#225EA8") +
        ggplot2::geom_vline(xintercept = reference_length + 0.5, linetype = "dashed") +
        ggplot2::annotate(
          "text",
          x = (reference_length + 1) / 2,
          y = Inf,
          label = "Phase I: reference",
          vjust = 1.4
        ) +
        ggplot2::annotate(
          "text",
          x = reference_length + (length(x$new_data) + 1) / 2,
          y = Inf,
          label = "Phase II: monitoring",
          vjust = 1.4
        ) +
        ggplot2::geom_point(
          data = plot_data[plot_data$signal, , drop = FALSE],
          color = "#D7301F",
          size = 2.5
        ) +
        ggplot2::labs(
          x = "Time",
          y = "Observed value",
          title = "QMM reference and monitoring observations",
          subtitle = "The dashed line separates Phase I and Phase II"
        ) +
        ggplot2::theme_minimal(base_size = 12)
    )
  }

  panel_levels <- c("Observed series", "QMM statistic")
  observed_data <- transform(
    plot_data,
    panel = factor("Observed series", levels = panel_levels)
  )
  statistic_data <- data.frame(
    time = reference_length + x$results$horizon,
    value = x$results$statistic,
    signal = x$results$signal,
    panel = factor("QMM statistic", levels = panel_levels)
  )
  limit_data <- data.frame(
    time = reference_length + x$results$horizon,
    value = x$results$control_limit,
    panel = factor("QMM statistic", levels = panel_levels)
  )
  overview_data <- rbind(observed_data, statistic_data)
  phase_labels <- data.frame(
    time = rep(
      c(
        (reference_length + 1) / 2,
        reference_length + (length(x$new_data) + 1) / 2
      ),
      times = length(panel_levels)
    ),
    value = Inf,
    label = c(
      "Phase I: reference",
      "Phase II: monitoring",
      "Phase I: model fitting",
      "Phase II: monitoring"
    ),
    panel = factor(rep(panel_levels, each = 2L), levels = panel_levels)
  )
  phase1_statistic_region <- data.frame(
    panel = factor("QMM statistic", levels = panel_levels)
  )

  ggplot2::ggplot(overview_data, ggplot2::aes(x = time, y = value)) +
    ggplot2::geom_rect(
      data = phase1_statistic_region,
      xmin = 0.5,
      xmax = reference_length + 0.5,
      ymin = -Inf,
      ymax = Inf,
      fill = "#F2F2F2",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_line(
      data = observed_data,
      color = "#225EA8"
    ) +
    ggplot2::geom_line(
      data = statistic_data,
      color = "#225EA8"
    ) +
    ggplot2::geom_line(
      data = limit_data,
      linetype = "dashed"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = signal),
      size = 2.2
    ) +
    ggplot2::geom_vline(
      xintercept = reference_length + 0.5,
      linetype = "dotted"
    ) +
    ggplot2::geom_text(
      data = phase_labels,
      vjust = 1.4,
      size = 3.5,
      inherit.aes = FALSE,
      mapping = ggplot2::aes(x = time, y = value, label = label)
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(panel),
      scales = "free_y"
    ) +
    ggplot2::scale_color_manual(
      values = c(`FALSE` = "#225EA8", `TRUE` = "#D7301F")
    ) +
    ggplot2::labs(
      x = "Time",
      y = NULL,
      color = "Signal",
      title = "QMM Phase I and Phase II overview",
      subtitle = "QMM statistics and pointwise limits are defined for Phase II"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}
