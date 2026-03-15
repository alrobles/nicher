# get_optim_par Get the parameters for the optimization of the Mahalanobis ellipse

get_optim_par Get the parameters for the optimization of the Mahalanobis
ellipse

## Usage

``` r
get_optim_par(df)
```

## Arguments

- df:

  A data frame with environmental information

## Value

A list with two objects. A vector with centers of ellipse and a matrix
with the inverse of covariance matrix

## Examples

``` r
get_optim_par(spOccPnts)
#> $mu
#>     bio1WH    bio12WH 
#>   21.39775 1872.49315 
#> 
#> $A
#>             bio1WH       bio12WH
#> [1,]  0.1312938652 -1.165306e-04
#> [2,] -0.0001165306  2.284606e-06
#> 
```
