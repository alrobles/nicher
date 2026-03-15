# Get ellipsoid parameters. A function to compute average and the inverse of covariance matrix from environmental data

Get ellipsoid parameters. A function to compute average and the inverse
of covariance matrix from environmental data

## Usage

``` r
get_ellip_par(env)
```

## Arguments

- env:

  A data frame containing environmental variables

## Value

A list with computed average of environmental variables and the
covariance matrix

## Examples

``` r
get_ellip_par(spOccPnts)
#> $mu
#>     bio1WH    bio12WH 
#>   21.39775 1872.49315 
#> 
#> $S
#>             bio1WH     bio12WH
#> bio1WH    7.977662    406.9154
#> bio12WH 406.915434 458467.6979
#> 
```
