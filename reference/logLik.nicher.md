# Log-likelihood of a fitted nicher model

Returns the \*\*un-penalised\*\* log-likelihood of a `"nicher"` fit at
the converged \\\theta\\. This is the appropriate quantity for
likelihood-based model comparison via
[`AIC`](https://rdrr.io/r/stats/AIC.html) or
[`BIC`](https://rdrr.io/r/stats/AIC.html): even when
[`optimize_niche`](https://alrobles.github.io/nicher/reference/optimize_niche.md)
minimises a ridge-penalised objective \\-\log L(\theta) +
\mathrm{penalty}(\theta)\\, the comparison metric is the bare \\\log
L(\theta)\\ evaluated at the optimised \\\theta\\.

## Usage

``` r
# S3 method for class 'nicher'
logLik(object, ...)
```

## Arguments

- object:

  A `"nicher"` object returned by
  [`optimize_niche`](https://alrobles.github.io/nicher/reference/optimize_niche.md).

- ...:

  Ignored.

## Value

An object of class `"logLik"` with attributes `df` (number of free
parameters) and `nobs` (number of occurrence rows used in the fit).

## Details

For penalised fits this means the value returned is \*\*not\*\* the same
as `x$best$loglik`: `best$loglik` is the negative of the penalised
objective at the optimum (and is monotone in fit quality but not
directly comparable across penalty strengths), while `logLik(x)` is the
bare data log-likelihood. `logLik(x)` is what AIC/BIC require.

Stored on `x$best$loglik_unpenalised` by
[`optimize_niche`](https://alrobles.github.io/nicher/reference/optimize_niche.md);
older fits without that field are detected and an informative error is
raised.

## Effective degrees of freedom

`df` is reported naively as `length(theta)`: 2p + p(p-1)/2 for the
Gaussian families, 3p + p(p-1)/2 for the skew-normal families, 3p +
p(p-1)/2 + 1 for the skew-t families. This \*\*ignores\*\* the shrinkage
induced by the ridge penalties (`prior_mu_lambda`,
`prior_log_sigma_lambda`, `prior_alpha_lambda`). For weak penalties
(\\\lambda \le 1\\) the bias is small; for stronger penalties the
effective df is lower than `length(theta)` so the penalty (\\2k\\ or \\k
\log n\\) is conservative — IC values will under-favour the
more-regularised model. Quantifying effective df via the trace of the
influence matrix is left to a future release.

## Examples

``` r
# \donttest{
fit <- optimize_niche(env_occ = example_env_occ_2d,
                      env_m   = example_env_m_2d,
                      num_starts = 5L, likelihood = "weighted")
logLik(fit)
#> 'log Lik.' -650.4004 (df=5)
AIC(fit)
#> [1] 1310.801
BIC(fit)
#> [1] 1322.253
# }
```
