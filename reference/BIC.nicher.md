# Bayesian information criterion for a nicher fit

Computes \\-2 \log L + k \log n\\, where \\n\\ is `nobs(object)`. Calls
[`BIC`](https://rdrr.io/r/stats/AIC.html).

## Usage

``` r
# S3 method for class 'nicher'
BIC(object, ...)
```

## Arguments

- object:

  a fitted model object for which there exists a `logLik` method to
  extract the corresponding log-likelihood, or an object inheriting from
  class `logLik`.

- ...:

  Additional fitted model objects.

## Value

Numeric scalar, or a `data.frame` when multiple objects are passed.

## Examples

``` r
# \donttest{
fit <- optimize_niche(env_occ = example_env_occ_2d,
                      env_m   = example_env_m_2d,
                      num_starts = 5L, likelihood = "weighted")
BIC(fit)
#> [1] 1322.253
# }
```
