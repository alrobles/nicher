# Number of observations used to fit a nicher model

Returns the number of presence rows (`nrow(env_occ)`) that were passed
to
[`optimize_niche`](https://alrobles.github.io/nicher/reference/optimize_niche.md).
This is the sample size used by
[`BIC`](https://rdrr.io/r/stats/AIC.html) and is the conventional choice
for the presence-only and weighted families: even though the weighted
likelihood includes a denominator over background points, the observed
quantities being modeled are the \\n\_{occ}\\ occurrences.

## Usage

``` r
# S3 method for class 'nicher'
nobs(object, ...)
```

## Arguments

- object:

  A `"nicher"` object.

- ...:

  Ignored.

## Value

Integer scalar.

## Examples

``` r
# \donttest{
fit <- optimize_niche(env_occ = example_env_occ_2d,
                      env_m   = example_env_m_2d,
                      num_starts = 5L, likelihood = "weighted")
nobs(fit)
#> [1] 73
# }
```
