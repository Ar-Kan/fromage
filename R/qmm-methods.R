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

#' Select horizons for the bootstrap density panels of a plot
#'
#' When `requested` is `NULL`, picks at most `max_default` horizons that have
#' already been monitored (no greater than `monitored_length`), evenly spread
#' across the eligible range. An explicit `requested` may name any horizon
#' present in the calibration, including horizons beyond the current
#' monitoring progress; it is only checked against the calibrated range.
#'
#' @param calibrated_horizons Horizons present in a `qmm_calibration` object's
#'   `limits` (i.e. `calibration$limits$horizon`).
#' @param monitored_length Number of Phase II observations monitored so far
#'   (i.e. `nrow(monitor$results)`).
#' @param requested Optional whole-number horizon(s) to use instead of the
#'   default selection.
#' @param max_default Maximum number of horizons picked when `requested` is
#'   `NULL`.
#'
#' @return A sorted, de-duplicated integer vector of horizons.
#' @keywords internal
select_density_horizons <- function(calibrated_horizons, monitored_length, requested = NULL, max_default = 4L) {
  calibrated_horizons <- sort(unique(as.integer(calibrated_horizons)))

  if (is.null(requested)) {
    eligible <- calibrated_horizons[calibrated_horizons <= monitored_length]
    if (length(eligible) == 0L) {
      stop(
        "No calibrated horizon has been monitored yet; supply density_horizons explicitly.",
        call. = FALSE
      )
    }
    if (length(eligible) <= max_default) {
      return(eligible)
    }
    positions <- unique(round(seq(1L, length(eligible), length.out = max_default)))
    return(sort(eligible[positions]))
  }

  requested <- vapply(
    requested,
    function(value) validate_whole_number(value, "density_horizons"),
    integer(1)
  )
  requested <- sort(unique(requested))
  missing <- setdiff(requested, calibrated_horizons)
  if (length(missing) > 0L) {
    stop(
      "density_horizons is not part of the calibration: ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  requested
}

#' Build the Phase I/Phase II annotated series panel shared by several plots
#'
#' @param plot_data A data frame with `time`, `value`, and `signal` columns.
#' @param reference_length Number of Phase I observations.
#' @param new_data_length Number of Phase II observations.
#'
#' @return A `ggplot` object without a title or subtitle.
#' @keywords internal
qmm_series_panel <- function(plot_data, reference_length, new_data_length) {
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
      x = reference_length + (new_data_length + 1) / 2,
      y = Inf,
      label = "Phase II: monitoring",
      vjust = 1.4
    ) +
    ggplot2::geom_point(
      data = plot_data[plot_data$signal, , drop = FALSE],
      color = "#D7301F",
      size = 2.5
    ) +
    ggplot2::labs(x = "Time", y = "Observed value") +
    ggplot2::theme_minimal(base_size = 12)
}

#' Build the Phase II QMM statistic panel shared by several plots
#'
#' @param results A monitoring result's `results` data frame.
#' @param show_limit_legend Whether to map the control-limit line to a
#'   "Control limit" legend entry instead of drawing it as a bare dashed line.
#'
#' @return A `ggplot` object without a title or subtitle.
#' @keywords internal
qmm_statistic_panel <- function(results, show_limit_legend = FALSE) {
  panel <- ggplot2::ggplot(results, ggplot2::aes(x = horizon))
  panel <- if (show_limit_legend) {
    panel + ggplot2::geom_line(ggplot2::aes(y = control_limit, linetype = "Control limit"))
  } else {
    panel + ggplot2::geom_line(ggplot2::aes(y = control_limit), linetype = "dashed")
  }

  panel +
    ggplot2::geom_line(ggplot2::aes(y = statistic), color = "#225EA8") +
    ggplot2::geom_point(ggplot2::aes(y = statistic, color = signal), size = 2.2) +
    ggplot2::scale_color_manual(values = c(`FALSE` = "#225EA8", `TRUE` = "#D7301F")) +
    (if (show_limit_legend) {
      ggplot2::scale_linetype_manual(name = NULL, values = c(`Control limit` = "dashed"))
    }) +
    ggplot2::labs(x = "Monitoring horizon", y = "QMM statistic", color = "Signal") +
    ggplot2::theme_minimal(base_size = 12)
}

#' Build one bootstrap-statistic density panel
#'
#' @param values Bootstrap statistics simulated for one horizon.
#' @param horizon The horizon these values were simulated for.
#' @param control_limit The horizon's control limit.
#' @param alpha The calibration's exceedance probability (`1 -
#'   confidence_level`).
#' @param x_range A length-two numeric vector shared across every density
#'   panel in the same row, so the panels stay visually comparable.
#' @param show_y_label Whether to label the y-axis `"Density"`.
#'
#' @return A `ggplot` object.
#' @keywords internal
qmm_density_panel <- function(values, horizon, control_limit, alpha, x_range, show_y_label) {
  ggplot2::ggplot(data.frame(statistic = values), ggplot2::aes(x = statistic)) +
    ggplot2::geom_density(fill = "#D7301F", color = "#D7301F", alpha = 0.15) +
    ggplot2::geom_vline(xintercept = control_limit, linetype = "dotted") +
    ggplot2::annotate(
      "text",
      x = control_limit,
      y = -Inf,
      label = sprintf("Limit at \u03b1 = %s", format(signif(alpha, 3), trim = TRUE)),
      angle = 90,
      hjust = 0,
      vjust = -0.6,
      size = 2.6
    ) +
    ggplot2::coord_cartesian(xlim = x_range) +
    ggplot2::labs(
      x = "QMM statistic",
      y = if (show_y_label) "Density" else NULL,
      title = sprintf("Horizon %d", horizon)
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 10),
      panel.background = ggplot2::element_rect(fill = "#F2F2F2", color = NA)
    )
}

#' Plot a multimodel Q monitoring result
#'
#' @param x A monitoring result created by [monitor_qmm()] or [qmm_chart()].
#' @param type Plot the Phase II QMM statistic and control limit
#'   (`"statistic"`), the complete Phase I/Phase II series (`"series"`), an
#'   aligned two-panel view containing both (`"overview"`), or the same
#'   two-panel view with a row of bootstrap-statistic density panels in place
#'   of the Phase I region (`"bootstrap"`). `"bootstrap"` requires that `x`'s
#'   calibration was created with `calibrate_qmm(..., keep_bootstrap = TRUE)`.
#' @param density_horizons For `type = "bootstrap"`, the horizon(s) to draw a
#'   density panel for. Defaults to `NULL`, which auto-selects at most 4
#'   horizons already reached by monitoring (see [select_density_horizons()]).
#'   An explicit value may also name a calibrated horizon beyond the current
#'   monitoring progress.
#' @param ... Unused.
#'
#' @return A `ggplot` object, or (`type = "bootstrap"`) a `patchwork` object.
#' @export
plot.qmm_monitor <- function(x,
                             type = c("statistic", "series", "overview", "bootstrap"),
                             density_horizons = NULL,
                             ...) {
  type <- match.arg(type)

  if (type == "statistic") {
    return(
      qmm_statistic_panel(x$results) +
        ggplot2::labs(
          title = "Multimodel Q (QMM) control chart",
          subtitle = "Phase II; the dashed line is the pointwise upper control limit"
        )
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
      qmm_series_panel(plot_data, reference_length, length(x$new_data)) +
        ggplot2::labs(
          title = "QMM reference and monitoring observations",
          subtitle = "The dashed line separates Phase I and Phase II"
        )
    )
  }

  if (type == "bootstrap") {
    calibration <- x$calibration
    if (is.null(calibration$bootstrap)) {
      stop(
        "type = \"bootstrap\" requires bootstrap draws; ",
        "recalibrate with calibrate_qmm(..., keep_bootstrap = TRUE).",
        call. = FALSE
      )
    }

    horizons <- select_density_horizons(
      calibration$limits$horizon,
      monitored_length = nrow(x$results),
      requested = density_horizons
    )
    horizon_columns <- match(horizons, calibration$limits$horizon)
    control_limits <- calibration$limits$control_limit[horizon_columns]
    alpha <- 1 - calibration$confidence_level
    bootstrap_values <- lapply(horizon_columns, function(column) calibration$bootstrap[, column])
    x_range <- range(unlist(bootstrap_values), control_limits)

    density_plots <- Map(
      qmm_density_panel,
      values = bootstrap_values,
      horizon = horizons,
      control_limit = control_limits,
      MoreArgs = list(alpha = alpha, x_range = x_range),
      show_y_label = c(TRUE, rep(FALSE, length(horizons) - 1L))
    )

    density_row <- patchwork::wrap_elements(
      patchwork::wrap_plots(density_plots, nrow = 1L)
    )
    top_panel <- qmm_series_panel(plot_data, reference_length, length(x$new_data))
    statistic_panel <- qmm_statistic_panel(x$results, show_limit_legend = TRUE)
    bottom_row <- density_row + statistic_panel +
      patchwork::plot_layout(widths = c(reference_length, length(x$new_data)))

    return(
      (top_panel / bottom_row) +
        patchwork::plot_layout(guides = "collect") +
        patchwork::plot_annotation(
          title = "QMM Phase I and Phase II overview",
          subtitle = "QMM statistics and pointwise limits are defined for Phase II"
        )
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
