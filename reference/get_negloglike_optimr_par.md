# get_negloglike_optim_par Function to optime negloglike function given parameters from presence points and a sample of environmental points from a M hypothesis. Internally runs Returns a list of parameters

get_negloglike_optim_par Function to optime negloglike function given
parameters from presence points and a sample of environmental points
from a M hypothesis. Internally runs Returns a list of parameters

## Usage

``` r
get_negloglike_optimr_par(
  env_pts,
  M_pts,
  lower = FALSE,
  itnmax = 100,
  fastmode = FALSE,
  return_best = TRUE,
  precision = 3
)
```

## Arguments

- env_pts:

  A dataframe with environmental variables extracted from presence
  points

- M_pts:

  A dataframe with a sample of environmental values inside an M
  hypothesis region.

- lower:

  logical if is set create a low bowndarie for parameters from M points

- itnmax:

  The maximum number of iterations.

- fastmode:

  logical if is true the L-BFGS-B method is by default

- return_best:

  logical. If true return the best optimization. If false return all the
  possible optimizations

- precision:

  numeric. The decimal precision of the output variables. The default is
  3

## Value

A list with optimized parameters for negloglike function

## Examples

``` r
get_negloglike_optimr_par(head(spOccPnts, 10), samMPts, fastmode = TRUE, itnmax = 1)
#> fastmode
#> TRUE
#> Optimization starting 
#> Methods used:
#> L-BFGS-B
#> Optimization finished. 
#> L-BFGS-B is the returned method
#> $mu
#>   bio1WH  bio12WH 
#>   20.971 2019.703 
#> 
#> $S
#>          bio1WH    bio12WH
#> bio1WH    2.193    252.165
#> bio12WH 252.165 230331.789
#> 
```
