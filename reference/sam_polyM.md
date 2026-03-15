# Sample environmental points from M hypothesis

Sample environmental points from M hypothesis

## Usage

``` r
sam_polyM(M.shp, env, N = 100)
```

## Arguments

- M.shp:

  A shapefile with an M hypothesis of the species

- env:

  an stack of raster with environmental variables to extract information

- N:

  a numeric with the number of points of the sample, default is 10 000

## Value

a data frame with a sample of environmental values inside the M polygon

## Examples

``` r
M_path <- system.file("extdata", "Mshp_test.rds", package="nicher")
Mshp <- terra::unwrap(readr::read_rds(M_path))
stack_path <- system.file("extdata", "stack_1_12_crop.rds", package="nicher")
stack_1_12 <- terra::unwrap(readr::read_rds(stack_path))
sam_polyM(M.shp = Mshp, N = 5, env = stack_1_12)
#>        bio1WH bio12WH
#> [1,] 27.97160    1358
#> [2,] 19.08194    1830
#> [3,] 28.58866    1562
#> [4,] 25.78077    2003
#> [5,] 23.40339    1560
```
