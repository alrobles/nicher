# k-fold (or leave-one-out) cross-validation for a nicher fit

Refits the same likelihood family + same penalty knobs as the supplied
`nicher` object on each train fold of `env_occ`, warm-started from the
full-data \\\hat\theta\\, then scores the held-out fold by evaluating
the un-penalised log-likelihood at the fold's converged
\\\hat\theta\_{train}\\. Returns the \*\*summed\*\* held-out
log-likelihood across folds, which is the standard out-of-sample
log-likelihood used to tune ridge regularisation.

## Usage

``` r
cv_nicher(
  fit,
  env_occ,
  env_m = NULL,
  type = c("kfold", "loo"),
  k = 5L,
  seed = NULL,
  num_starts_cv = 1L,
  verbose = FALSE,
  ucminf_control = list(maxeval = 500L)
)
```

## Arguments

- fit:

  A `"nicher"` object returned by
  [`optimize_niche`](https://alrobles.github.io/nicher/reference/optimize_niche.md).

- env_occ:

  A `matrix` or `data.frame` of occurrences, identical (up to row order)
  to the one passed to the original
  [`optimize_niche()`](https://alrobles.github.io/nicher/reference/optimize_niche.md)
  call. Validated against `fit$best$env_occ_fingerprint`.

- env_m:

  A `matrix` / `data.frame` of background points, or `NULL` for
  presence-only families. Validated against
  `fit$best$env_m_fingerprint`.

- type:

  Character. `"kfold"` (default) or `"loo"`.

- k:

  Integer. Number of folds when `type = "kfold"`; ignored for `"loo"`.
  Default `5L`.

- seed:

  Integer or `NULL`. Random seed for fold assignment (`type = "kfold"`
  only). Default `NULL` (no seed).

- num_starts_cv:

  Integer. Number of starts per fold. Default `1L` (warm-start from
  `fit$best$theta`). Set higher for robustness in `"kfold"` mode; for
  `"loo"` this should essentially always be `1L` since you are running
  n_occ refits.

- verbose:

  Logical. Print per-fold progress.

- ucminf_control:

  Optional list passed to
  [`ucminf::ucminf`](https://rdrr.io/pkg/ucminf/man/ucminf.html) for
  each fold's optimiser; defaults to `list(maxeval = 500L)`.

## Value

A list with class `"nicher_cv"`:

- `cv_loglik`:

  Summed held-out log-likelihood across folds.

- `cv_loglik_mean`:

  Per-occurrence average: `cv_loglik / nobs`.

- `per_fold`:

  `data.frame(fold, n_test, loglik_test, convergence)`.

- `type`:

  `"kfold"` or `"loo"`.

- `k`:

  Number of folds (`n_occ` for LOO).

- `seed`:

  Seed used (or `NA_integer_`).

## Details

For penalised fits the un-penalised log-likelihood at the converged
\\\hat\theta\\ (what
[`logLik.nicher`](https://alrobles.github.io/nicher/reference/logLik.nicher.md)
returns) is monotone in the penalty strengths – the in-sample loglik is
always best at \\\lambda = 0\\. Cross-validation breaks that
monotonicity by scoring on data the optimiser did not see, so a properly
tuned penalty achieves higher out-of-sample loglik than the unpenalised
MLE.

## Caveats

\* Random k-fold CV assumes the occurrences are exchangeable. For
spatially autocorrelated data this can be optimistic; spatial-block CV
is left to a future release. \* LOO with `n_occ` large (e.g. Vicugna's
545) is feasible but non-trivial: 545 ucminf refits at ~0.1-1 s each =
~1-10 minutes. Consider `type = "kfold", k = 10L` as a faster proxy. \*
Each fold's refit uses
[`ucminf::ucminf`](https://rdrr.io/pkg/ucminf/man/ucminf.html) starting
from `fit$best$theta`. If the per-fold likelihood is not unimodal (rare
on dense data, common with skew families on sparse data), bump
`num_starts_cv` or use the original
[`optimize_niche()`](https://alrobles.github.io/nicher/reference/optimize_niche.md)
multistart machinery.

## See also

[`compare_nicher`](https://alrobles.github.io/nicher/reference/compare_nicher.md)
for AIC/BIC-based comparisons.

## Examples

``` r
# \donttest{
fit_w <- optimize_niche(
  env_occ = example_env_occ_2d, env_m = example_env_m_2d,
  num_starts = 5L, breadth = 0.45, likelihood = "weighted",
  prior_mu_lambda = 1.0
)
cv_nicher(fit_w, example_env_occ_2d, example_env_m_2d,
          type = "kfold", k = 5L, seed = 42L)
#> -- nicher_cv result --
#>   Type     : kfold (k = 5) 
#>   Family   : weighted 
#>   Folds    : 5 
#>   CV loglik: -658.526 (per-occ: -9.0209 )
# }
```
