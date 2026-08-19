# fromage

`fromage` implements the multimodel Q control chart, written as $Q_{MM}$, for monitoring changes in the serial
dependence of an autocorrelated process.

## Installation

Install the development version from GitHub with:

```r
# install.packages("pak")
pak::pak("Ar-Kan/fromage")
```

## Minimal usage

`fromage` implements the multimodel Q control chart, written as $Q_{MM}$, for monitoring changes in the serial
dependence of an autocorrelated process.

`process_series` must be one numeric vector or univariate time series ordered by collection time. In this example, its
first 100 observations are the in-control reference period.

```r
library(fromage)

chart <- qmm_chart(
  x = process_series,
  reference_length = 100
)

chart
signals(chart)
plot(chart, type = "overview")
```

## Documentation

- `vignette("qmm-chart", package = "fromage")` explains the methodology and follows a realistic simulated process from
  Phase I fitting through sequential monitoring.
- `vignette("qmm-simulation-study", package = "fromage")` reproduces the QMM portion of the original simulation-study
  design using the package API.
- `vignette("qmm-internals", package = "fromage")` documents the public and internal implementation functions.
- [`CONTRIBUTING.md`](CONTRIBUTING.md) describes the repository architecture, object contracts, tests, and development
  workflow.

## References

- Ljung, G. M. and Box, G. E. P. (1978). On a measure of lack of fit in time series models. *Biometrika*, 65 (2),
  297–303.
  <https://doi.org/10.1093/biomet/65.2.297>
- Burnham, K. P. and Anderson, D. R. (2002). *Model Selection and Multimodel Inference: A Practical
  Information-Theoretic Approach*. Springer.
