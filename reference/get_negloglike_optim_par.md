# get_negloglike_optim_par Function to optime negloglike function given parameters from presence points and a sample of environmental points from a M hypothesis. Returns a list of

get_negloglike_optim_par Function to optime negloglike function given
parameters from presence points and a sample of environmental points
from a M hypothesis. Returns a list of

## Usage

``` r
get_negloglike_optim_par(env_pts, M_pts)
```

## Arguments

- env_pts:

  A dataframe with environmental variables extracted from presence
  points

- M_pts:

  A dataframe with a sample of environmental values inside an M
  hypothesis region.

## Value

A list with optimized parameters for negloglike function. mu is a vector
of centroids and S is a covariance matrix.

## Examples

``` r
get_negloglike_optim_par(head(spOccPnts, 30), samMPts)
#> $mu
#>     bio1WH    bio12WH 
#>   17.17463 2414.46048 
#> 
#> $S
#>             bio1WH    bio12WH
#> bio1WH     8.57205  -2314.157
#> bio12WH 2009.52866 465713.414
#> 
```
