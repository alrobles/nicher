# Negative log likelyhood

Negative log likelyhood

## Usage

``` r
negloglike_multivariable(mu, S, sam1, sam2)
```

## Arguments

- mu:

  A vector mu of parameters

- S:

  The covariance matrix from environmental data frame

- sam1:

  A data.frame containing the original sample of environmental
  combinations that correspond to presences

- sam2:

  A data.frame containing a second random sample of environmental
  combinations which come from the area of study (M)

## Value

A negative log likelihood value

## Examples

``` r
par <- get_ellip_par(spOccPnts)
negloglike_multivariable(par$mu, par$S, spOccPnts, samMPts)
#> [1] 661.8376
```
