# Calculates the inverse of covariance matrix from environmental variables

Calculates the inverse of covariance matrix from environmental variables

## Usage

``` r
get_A_matrix(df)
```

## Arguments

- df:

  A matrix or data frame with environmental variables

## Value

A covariance matrix

## Examples

``` r
get_A_matrix(spOccPnts)
#>             bio1WH       bio12WH
#> [1,]  0.1312938652 -1.165306e-04
#> [2,] -0.0001165306  2.284606e-06
```
