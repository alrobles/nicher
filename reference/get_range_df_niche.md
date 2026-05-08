# Generate reasonable ranges for starting parameters (math scale)

Ensures that env_data is strictly numeric before computing ranges. This
prevents downstream issues in Sobol designs and multi-start optimization
(e.g., character / factor leakage).

## Usage

``` r
get_range_df_niche(
  env_data,
  quant_vec = c(0.1, 0.5, 0.9),
  skew = FALSE,
  skew_t = FALSE
)
```

## Arguments

- env_data:

  Numeric matrix or data frame of environmental values.

- quant_vec:

  Numeric vector of quantiles (length 3).

- skew:

  Logical. If `TRUE`, append `p` skew parameters
  (`alpha_1, ..., alpha_p`) to the parameter vector with the default
  Sobol range `[-3, 0, 3]`. Default `FALSE` keeps the Gaussian-only
  parameterization.

- skew_t:

  Logical. If `TRUE`, append a single `log_r` parameter
  (degrees-of-freedom for the multivariate non-central skew-t) AFTER the
  alpha block. The Sobol range is `[log 2, log 10, log 100]`; the
  optimizer is then free to move `log_r` anywhere on the real line.
  Implies `skew = TRUE`.

## Value

A data frame with rows = parameter names and columns lower, center,
upper.
