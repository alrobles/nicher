# Akaike information criterion for a nicher fit

Computes \\-2 \log L + 2 k\\, where \\\log L\\ is the un-penalised
log-likelihood (see
[`logLik.nicher`](https://alrobles.github.io/nicher/reference/logLik.nicher.md))
and \\k\\ is the number of free parameters.

## Usage

``` r
# S3 method for class 'nicher'
AIC(object, ..., k = 2)
```

## Arguments

- object:

  a fitted model object for which there exists a `logLik` method to
  extract the corresponding log-likelihood, or an object inheriting from
  class `logLik`.

- ...:

  Additional fitted model objects.

- k:

  Numeric, the penalty per parameter; default `2` (AIC). Pass
  `k = log(nobs(fit))` to recover BIC.

## Value

Numeric scalar, or a `data.frame` when multiple objects are passed.

## Examples

``` r
# \donttest{
fit <- optimize_niche(env_occ = example_env_occ_2d,
                      env_m   = example_env_m_2d,
                      num_starts = 5L, likelihood = "weighted")
AIC(fit)
#> [1] 1310.801
# }
```
