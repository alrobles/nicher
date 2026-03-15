# Ellipse constructor

Constructor function that create an ellipse.

## Usage

``` r
ellipse(param_list)
```

## Arguments

- param_list:

  A list of parameters of an ellipse object. First element is the mu
  vector. Second element is the square matrix of covariates

## Value

`ellipse`

## Examples

``` r
ellipse(get_ellip_par(spOccPnts))
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
