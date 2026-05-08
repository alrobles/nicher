# Cross-Validation: Penalty Tuning and Model Selection

## Overview

Ridge penalties stabilise niche estimates when sample sizes are small or
environmental variables are correlated. `nicher` supports three penalty
knobs:

| Parameter                | Controls                           | Default |
|--------------------------|------------------------------------|---------|
| `prior_log_sigma_lambda` | Shrinkage on $`\log\sigma`$        | `0`     |
| `prior_mu_lambda`        | Shrinkage on $`\mu`$ toward center | `0`     |
| `prior_alpha_lambda`     | Shrinkage on skewness $`\alpha`$   | `0`     |

Setting a penalty too high over-regularises (biased niche); setting it
too low under-regularises (high variance). Cross-validation selects the
penalty strength that maximises out-of-sample log-likelihood.

## 1 Fit a baseline model

``` r

library(nicher)
set.seed(42)

fit_base <- optimize_niche(
  env_occ    = example_env_occ_2d,
  env_m      = example_env_m_2d,
  likelihood = "weighted",
  num_starts = 5L,
  breadth    = 0.45,
  prior_mu_lambda = 0,
  verbose    = FALSE
)
fit_base
#> -- nicher optimization result --
#>   Likelihood : weighted 
#>   Starts     : 6 (converged: 6 )
#>   Best loglik: -645.0662 
#>   Convergence: 1 
#>   eta        : 1 
#>   Variables  : bio1WH, bio12WH
```

## 2 Run k-fold cross-validation

[`cv_nicher()`](https://alrobles.github.io/nicher/reference/cv_nicher.md)
refits the same model on each train fold (warm-started from the
full-data optimum) and scores the held-out fold using the un-penalised
log-likelihood.

``` r

cv_base <- cv_nicher(
  fit     = fit_base,
  env_occ = example_env_occ_2d,
  env_m   = example_env_m_2d,
  type    = "kfold",
  k       = 5L,
  seed    = 42L,
  verbose = FALSE
)
cv_base
#> -- nicher_cv result --
#>   Type     : kfold (k = 5) 
#>   Family   : weighted 
#>   Folds    : 5 
#>   CV loglik: -650.923 (per-occ: -8.9168 )
```

## 3 Compare penalty strengths

Fit models at different penalty levels and compare their CV
log-likelihoods. Higher (less negative) CV log-likelihood indicates
better out-of-sample fit.

``` r

lambdas <- c(0, 0.5, 1.0, 2.0, 5.0)
cv_results <- data.frame(
  lambda   = lambdas,
  cv_loglik     = NA_real_,
  cv_loglik_mean = NA_real_
)

for (i in seq_along(lambdas)) {
  fit_i <- optimize_niche(
    env_occ    = example_env_occ_2d,
    env_m      = example_env_m_2d,
    likelihood = "weighted",
    num_starts = 5L,
    breadth    = 0.45,
    prior_mu_lambda = lambdas[i],
    verbose    = FALSE
  )
  cv_i <- cv_nicher(
    fit     = fit_i,
    env_occ = example_env_occ_2d,
    env_m   = example_env_m_2d,
    type    = "kfold",
    k       = 5L,
    seed    = 42L,
    verbose = FALSE
  )
  cv_results$cv_loglik[i]      <- cv_i$cv_loglik
  cv_results$cv_loglik_mean[i] <- cv_i$cv_loglik_mean
}

cv_results
#>   lambda cv_loglik cv_loglik_mean
#> 1    0.0 -650.9230      -8.916754
#> 2    0.5 -656.6578      -8.995313
#> 3    1.0 -658.5260      -9.020904
#> 4    2.0 -660.9557      -9.054188
#> 5    5.0 -665.4376      -9.115583
```

The penalty value that maximises `cv_loglik` is the best choice:

``` r

best_idx <- which.max(cv_results$cv_loglik)
cat("Best lambda:", cv_results$lambda[best_idx], "\n")
#> Best lambda: 0
cat("CV loglik  :", cv_results$cv_loglik[best_idx], "\n")
#> CV loglik  : -650.923
```

## 4 Cross-validation across likelihood families

CV can also compare different likelihood families. Refit each family
with the same data and compare their out-of-sample log-likelihoods:

``` r

families <- c("presence_only", "weighted", "ip_weighted")
cv_family <- data.frame(
  family         = families,
  cv_loglik      = NA_real_,
  cv_loglik_mean = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(families)) {
  fam <- families[i]
  needs_m <- fam %in% c("weighted", "ip_weighted")
  fit_i <- optimize_niche(
    env_occ    = example_env_occ_2d,
    env_m      = if (needs_m) example_env_m_2d else NULL,
    likelihood = fam,
    num_starts = 5L,
    breadth    = 0.45,
    verbose    = FALSE
  )
  cv_i <- cv_nicher(
    fit     = fit_i,
    env_occ = example_env_occ_2d,
    env_m   = if (needs_m) example_env_m_2d else NULL,
    type    = "kfold",
    k       = 5L,
    seed    = 42L,
    verbose = FALSE
  )
  cv_family$cv_loglik[i]      <- cv_i$cv_loglik
  cv_family$cv_loglik_mean[i] <- cv_i$cv_loglik_mean
}

cv_family
#>          family cv_loglik cv_loglik_mean
#> 1 presence_only -630.3632      -8.635113
#> 2      weighted -658.5260      -9.020904
#> 3   ip_weighted -640.1458      -8.769120
```

## 5 Leave-one-out cross-validation

For small datasets, leave-one-out (LOO) CV uses every observation as a
test fold exactly once. This is more expensive but has lower variance
than k-fold:

``` r

cv_loo <- cv_nicher(
  fit     = fit_base,
  env_occ = example_env_occ_2d,
  env_m   = example_env_m_2d,
  type    = "loo",
  verbose = FALSE
)
cv_loo
#> -- nicher_cv result --
#>   Type     : loo  
#>   Family   : weighted 
#>   Folds    : 73 
#>   CV loglik: -649.3506 (per-occ: -8.8952 )
```

## 6 Per-fold diagnostics

The `per_fold` data frame in the result object shows convergence status
and log-likelihood for each fold, useful for detecting problematic
subsets:

``` r

head(cv_base$per_fold)
#>   fold n_test loglik_test convergence
#> 1    1     15   -134.8284           1
#> 2    2     14   -133.4668           1
#> 3    3     15   -130.5017           4
#> 4    4     14   -126.9579           4
#> 5    5     15   -125.1682           4
```

Folds with `convergence != 1` may need more optimizer iterations
(increase `ucminf_control$maxeval`).

## Caveats

- **Spatial autocorrelation**: Random k-fold CV assumes exchangeable
  observations. For spatially autocorrelated data, the CV log-likelihood
  may be optimistic. Spatial-block CV is planned for a future release.
- **Computational cost**: LOO with large $`n`$ (e.g. 500+) runs one
  refit per observation. Consider `type = "kfold", k = 10L` as a faster
  proxy.
- **Different normalisation**: Presence-only and weighted families have
  different likelihood normalising constants, so their CV
  log-likelihoods are not strictly comparable. Use this comparison as a
  diagnostic rather than a formal test.

## References

Jimenez, L., Soberon, J., Christen, J. A., & Soto, D. (2022). Dealing
with critical uncertainty in ecological niche models. *Ecological
Modelling*, 468, 109823.
<https://doi.org/10.1016/j.ecolmodel.2021.109823>
