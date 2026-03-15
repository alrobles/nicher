# Get the negative log-likelihood value from quadratic forms of environmental information from presence points and a sample from an M hypothesis

Get the negative log-likelihood value from quadratic forms of
environmental information from presence points and a sample from an M
hypothesis

## Usage

``` r
get_negative_log(q1, q2)
```

## Arguments

- q1:

  quadratic terms of presence points

- q2:

  quadratic terms of M points

## Value

A single value of the negative log-likelihood

## Examples

``` r
par <- get_optim_par(spOccPnts)
q1 <- mahalanobis(spOccPnts, par$mu, par$A, inverted = TRUE)
q2 <- mahalanobis(samMPts, par$mu, par$A, inverted = TRUE)
get_negative_log(q1, q2)
#> [1] 661.8376
```
