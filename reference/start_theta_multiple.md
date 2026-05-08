# Generate multiple starting points for niche model optimization

Creates a set of starting parameter vectors on the math scale for use
with the log-likelihood functions. Ensures all starting values are
strictly numeric and finite, even if env_data is supplied as a data
frame.

## Usage

``` r
start_theta_multiple(
  env_data,
  num_starts = 100,
  quant_vec = c(0.1, 0.5, 0.9),
  method = "sobol",
  skew = FALSE,
  skew_t = FALSE
)
```

## Arguments

- env_data:

  Environmental data (matrix or data frame). Must be numeric.

- num_starts:

  Integer, number of starting points.

- quant_vec:

  Quantiles for mu ranges.

- method:

  "sobol" (Sobol design) or "uniform".

- skew:

  Logical. If `TRUE`, append `p` skewness parameters
  (`alpha_1, ..., alpha_p`) to each starting vector with the default
  range `[-3, 3]`. Use this for the `"skew_normal"` /
  `"skew_normal_weighted"` likelihoods. Default `FALSE`.

- skew_t:

  Logical. If `TRUE`, additionally append a single `log_r` parameter
  (default range `[log 2, log 100]`) after the alpha block. Use for
  `"skew_t"` / `"skew_t_weighted"`. Implies `skew = TRUE`.

## Value

A data frame of dimension num_starts × num_parameters.

## Examples

``` r
starts <- start_theta_multiple(
  env_data   = example_env_occ_2d,
  num_starts = 5L,
  method     = "uniform"
)
dim(starts)
#> [1] 5 5
```
