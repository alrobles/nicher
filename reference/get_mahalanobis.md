# Mahalanobis distance

Mahalanobis distance

## Usage

``` r
get_mahalanobis(df, el_pars)
```

## Arguments

- df:

  A data frame with environmental information according to presence
  points

- el_pars:

  A list with ellipsoid parameters

## Value

A vector with Mahalanobis distance according to each point in
environmental space

## Examples

``` r
pars <- get_negloglike_optimr_par(spOccPnts, samMPts, lower = TRUE, fastmode = TRUE, itnmax = 1 )
#> fastmode
#> TRUE
#> Optimization starting 
#> Methods used:
#> L-BFGS-B
#> Optimization finished. 
#> L-BFGS-B is the returned method
get_mahalanobis(spOccPnts, pars)
#>  [1] 0.07287030 0.99869159 0.09564457 0.63016263 2.01990473 0.51931505
#>  [7] 2.20077587 0.21490163 1.10913071 0.04071999 4.67876906 2.47257954
#> [13] 0.68555273 0.11949621 1.75405465 0.18425300 2.48942359 3.20842080
#> [19] 4.73627253 1.94276915 3.33718185 5.13770990 2.61189057 5.92324522
#> [25] 6.12884817 5.41393920 5.13991846 7.07393541 7.36947051 2.04047284
#> [31] 4.30919824 6.19868070 0.32145328 0.97714704 2.21604333 0.87478588
#> [37] 0.48012289 1.41254738 4.49472261 4.31987731 4.20598738 1.88259745
#> [43] 2.26235410 4.18590750 3.50784068 0.14877050 3.83386404 1.32431089
#> [49] 2.18676440 0.89803496 0.49888836 1.89048191 1.00330190 1.73112286
#> [55] 0.90010394 2.89855481 1.24481934 0.87478588 0.25573058 0.75401744
#> [61] 0.96569267 0.32145328 0.86791093 2.08410671 0.18425300 4.73339288
#> [67] 1.45887627 1.08380417 2.29960504 3.27393291 1.27317703 0.37206971
#> [73] 0.03582950
```
