# Introduction to nicher

## Overview

The `nicher` package estimates ecological niche models using ellipsoidal
geometry under an M hypothesis (Jimenez et al. 2022). It fits seven
likelihood families via multi-start optimisation with Sobol sequences,
backed by a compiled C++/Eigen engine for speed.

## Example data

The package ships with example datasets for *Vicugna vicugna*:

``` r

library(nicher)

str(example_env_occ_2d)
#>  num [1:73, 1:2] 20.7 22.4 22.3 19.3 23.9 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : chr [1:2] "bio1WH" "bio12WH"
str(example_env_m_2d)
#>  num [1:10000, 1:2] 25.5 26.7 21.5 22.1 22.4 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : chr [1:2] "bio1WH" "bio12WH"
```

## Fitting a niche model

The main entry point is
[`optimize_niche()`](https://alrobles.github.io/nicher/reference/optimize_niche.md).
Here we fit a presence-only model with a small number of starts for
speed:

``` r

set.seed(42)
fit_po <- optimize_niche(
  env_occ    = example_env_occ_2d,
  env_m      = NULL,
  num_starts = 5L,
  breadth    = 0.1,
  likelihood = "presence_only"
)
print(fit_po)
#> -- nicher optimization result --
#>   Likelihood : presence_only 
#>   Starts     : 5 (converged: 5 )
#>   Best loglik: -621.9007 
#>   Convergence: 1 
#>   eta        : 1 
#>   Variables  : bio1WH, bio12WH
```

## Convergence diagnostics

The [`assess()`](https://alrobles.github.io/nicher/reference/assess.md)
generic checks whether multiple starts converged to the same optimum:

``` r

diag <- assess(fit_po)
diag$flag
#> [1] "accepted_global"
diag$recommendation
#> [1] "Global optimum likely found: multiple converged starts agree closely."
```

## Extracting ellipsoid parameters

Recover the estimated mean vector and covariance matrix on the natural
scale:

``` r

pars <- get_ellipsoid_pars(example_env_occ_2d)
pars$mu
#>     bio1WH    bio12WH 
#>   21.39775 1872.49315
pars$s_mat
#>             bio1WH     bio12WH
#> bio1WH    7.977662    406.9154
#> bio12WH 406.915434 458467.6979
```

## Weighted model

For models that correct for sampling bias in the accessible area (M),
supply `env_m`:

``` r

set.seed(42)
fit_w <- optimize_niche(
  env_occ    = example_env_occ_2d,
  env_m      = example_env_m_2d,
  num_starts = 5L,
  breadth    = 0.1,
  likelihood = "weighted"
)
print(fit_w)
#> -- nicher optimization result --
#>   Likelihood : weighted 
#>   Starts     : 6 (converged: 6 )
#>   Best loglik: -654.8488 
#>   Convergence: 1 
#>   eta        : 1 
#>   Variables  : bio1WH, bio12WH
```

## Model comparison with AIC/BIC

Compare competing models using information criteria:

``` r

cmp <- compare_nicher(
  presence_only = fit_po,
  weighted      = fit_w
)
print(cmp[, c("model", "likelihood", "AIC", "dAIC", "BIC", "dBIC")])
#>           model    likelihood      AIC     dAIC      BIC     dBIC
#> 1 presence_only presence_only 1253.801  0.00000 1265.254  0.00000
#> 2      weighted      weighted 1310.801 56.99943 1322.253 56.99943
```

## References

Jimenez, L., Soberon, J., Christen, J. A., & Soto, D. (2022). Dealing
with critical uncertainty in ecological niche models. *Ecological
Modelling*, 468, 109823.
<https://doi.org/10.1016/j.ecolmodel.2021.109823>
