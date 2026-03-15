# get_ENM_par

get ellipsoid parameters for a Ecollogical fundamental niche model.

## Usage

``` r
get_ENM_par(
  occPts,
  env,
  M_shp = NULL,
  method = c("unbound", "bound", "mahalanobis"),
  fastmode = TRUE,
  itnmax = 100,
  return_best = TRUE
)
```

## Arguments

- occPts:

  A data frame with ocurrence points. Has species name, longitude and
  latitude

- env:

  A raster stack

- M_shp:

  A shape file as M hypothesis

- method:

  Method selection. By default, optimization of the negloglike function
  with lower bound

- fastmode:

  logical if is true the L-BFGS-B method is by default

- itnmax:

  The maximum number of iterations.

- return_best:

  logical. If true return the best optimization. If false return all the
  possible optimizations

## Value

a list of parameters

## Examples

``` r
M_path <- system.file("extdata", "Mshp_test.rds", package="nicher")
Mshp <- terra::unwrap(readr::read_rds(M_path))
stack_path <- system.file("extdata", "stack_1_12_crop.rds", package="nicher")
stack_1_12 <- terra::unwrap(readr::read_rds(stack_path))
get_ENM_par(rawSpOccPnts, stack_1_12, Mshp, method = "mahalanobis")
#> $mu
#>     bio1WH    bio12WH 
#>   21.41451 1882.60274 
#> 
#> $S
#>             bio1WH     bio12WH
#> bio1WH    8.073574    433.5456
#> bio12WH 433.545626 455709.4650
#> 
```
