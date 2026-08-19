# Developing fromage

This guide is intended for developers working from the repository. Users who want to apply the QMM chart should start
with `README.md` and
`vignette("qmm-chart", package = "fromage")`.

## Architectural overview

The public workflow has three explicit stages and one convenience wrapper:

```text
Phase I reference series
        |
        v
     fit_qmm() ----------> qmm_fit
                                |
                                v
                         calibrate_qmm() ----> qmm_calibration
                                                   |
Phase II observations -----------------------------+
                                                   v
                                            monitor_qmm()
                                                   |
                                                   v
                                              qmm_monitor
```

`qmm_chart()` performs all three stages when one vector already contains the Phase I reference observations followed by
the Phase II observations. Keeping the stages separate is useful when model fits and bootstrap diagnostics need to be
inspected or when observations arrive sequentially.

The statistical responsibilities are deliberately separated:

- `fit_qmm()` estimates candidate ARMA models and their information-criterion weights from Phase I only.
- `qmm_statistic()` filters a supplied window with the fixed Phase I models, calculates model-specific Ljung-Box
  statistics, and combines them using the Phase I weights.
- `calibrate_qmm()` estimates pointwise, horizon-specific upper limits with a parametric bootstrap.
- `monitor_qmm()` constructs the coupled window at each Phase II horizon and compares its statistic with the
  corresponding calibrated limit.

Do not move information from Phase II into model fitting or calibration unless the methodology itself is intentionally
being changed and revalidated.

## Source-file map

| Path                  | Responsibility                                                                                              |
|-----------------------|-------------------------------------------------------------------------------------------------------------|
| `R/fromage-package.R` | Package-level roxygen documentation.                                                                        |
| `R/validation.R`      | Input validation, local random seeds, warning capture, and ARMA admissibility checks.                       |
| `R/qmm-candidates.R`  | Candidate ARMA order construction and validation.                                                           |
| `R/qmm-fit.R`         | Phase I model fitting, information criteria, model inclusion, and weights.                                  |
| `R/qmm-statistic.R`   | Fixed-coefficient residuals and model-specific and multimodel Q statistics.                                 |
| `R/qmm-bootstrap.R`   | Coefficient draws, conditional ARMA simulation, single bootstrap replicates, and control-limit calibration. |
| `R/qmm-monitor.R`     | Batch and sequential monitoring, `qmm_chart()`, `signals()`, and `update()`.                                |
| `R/qmm-methods.R`     | Print, summary, and plot methods for QMM objects.                                                           |
| `R/arma-simulation.R` | Reproducible coefficient-change simulations used by examples and tests.                                     |
| `R/globals.R`         | Declarations needed to avoid false global-variable notes from plotting code.                                |

## Object contracts

The package uses small S3 objects. Functions that consume an object should use its class validator before relying on its
contents.

## Development workflow

Run commands from the package root. During interactive development:

```r
pkgload::load_all(".")
devtools::document(".")
devtools::test(".")
```

The repository scripts provide repeatable workflows at `scripts/`.
