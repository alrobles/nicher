# Compare nicher fits side-by-side via AIC and BIC

Builds a one-row-per-model summary table for a list of `"nicher"` fits,
including log-likelihood, parameter count, AIC, BIC, deltas relative to
the best model, and Akaike weights.

## Usage

``` r
compare_nicher(
  ...,
  sort_by = c("AIC", "BIC", "loglik"),
  comparison_basis = c("unpenalised", "penalised")
)
```

## Arguments

- ...:

  Two or more `"nicher"` objects, optionally named, OR a single named
  list of such objects.

- sort_by:

  Character. One of `"AIC"` (default), `"BIC"`, or `"loglik"`: which
  metric to order rows by. Best model is always at the top.

- comparison_basis:

  Character, one of `"unpenalised"` (default) or `"penalised"`.
  Determines which log-likelihood is used as the basis for AIC and BIC.

  - `"unpenalised"`: the bare data log-likelihood \\\ell(\hat\theta)\\
    returned by
    [`logLik.nicher`](https://alrobles.github.io/nicher/reference/logLik.nicher.md).
    Standard frequentist convention, and the right basis when comparing
    \*\*different likelihood families\*\* at fixed penalty knobs.

  - `"penalised"`: the optimiser's actual objective,
    \\\ell(\hat\theta) - \mathrm{pen}(\hat\theta)\\, i.e.
    `x$best$loglik`. The right basis when comparing \*\*the same family
    at different penalty strengths\*\*, since the unpenalised
    log-likelihood is monotone in \\(\lambda\_{\mu},
    \lambda\_{\log\sigma}, \lambda\_{\alpha})\\ and would always favour
    \\\lambda = 0\\. Note: this score is NOT a true frequentist IC – it
    is the penalised-objective AIC, an effective-fit score.

## Value

A `data.frame` with one row per fit and columns:

- `model`:

  Name (or call index) of the fit.

- `likelihood`:

  Likelihood family (`x$likelihood`).

- `loglik`:

  Un-penalised log-likelihood at convergence.

- `df`:

  Number of free parameters.

- `nobs`:

  Sample size (occurrences).

- `AIC`, `BIC`:

  Information criteria.

- `dAIC`, `dBIC`:

  Deltas vs. the minimum.

- `weight_AIC`:

  Akaike weights, \\\exp(-\Delta_i / 2) / \sum_j \exp(-\Delta_j / 2)\\.

- `convergence`:

  ucminf convergence code.

## Details

All fits MUST have been produced from the same `env_occ`; any mismatch
in `nobs` or in the variable-name ordering raises an error, since IC
values are not comparable across heterogeneous data.

## Caveats

AIC and BIC compare \*likelihoods\* on the SAME observed data. When
comparing weighted vs. presence-only families, both families evaluate
\\\log L\\ on \\n\_{occ}\\ occurrences; the weighted family additionally
normalises by an `env_m` integral, so its likelihood is on a different
scale. In this implementation, `logLik(weighted)` returns the full
weighted-likelihood value (occurrence sum minus log-denominator), which
is the right thing for \*\*within-family\*\* comparison (e.g.,
`weighted` vs. `skew_normal_weighted`) but should be interpreted with
care when contrasting weighted with presence-only.

## Examples

``` r
# \donttest{
fit_w   <- optimize_niche(example_env_occ_2d, example_env_m_2d,
                          num_starts = 5L, likelihood = "weighted")
fit_skn <- optimize_niche(example_env_occ_2d, example_env_m_2d,
                          num_starts = 5L, likelihood = "skew_normal_weighted")
compare_nicher(weighted = fit_w, skew_normal = fit_skn)
#>         model           likelihood    loglik df nobs      AIC     dAIC      BIC
#> 1    weighted             weighted -650.4004  5   73 1310.801 0.000000 1322.253
#> 2 skew_normal skew_normal_weighted -649.0369  7   73 1312.074 1.272975 1328.107
#>       dBIC weight_AIC convergence
#> 1 0.000000  0.6539591           1
#> 2 5.853894  0.3460409           1
# }
```
