# Samples of environmental data from M hypothesis to estimate negative log likelihood from Vicugna vicugna.

Samples of environmental data from M hypothesis to estimate negative log
likelihood from Vicugna vicugna.

## Usage

``` r
example_vicugna
```

## Format

A list with 4 elements:

- env_m:

  Data frame with two columns and 10 000 rows each representing bio1 and
  bio12 of the accessibility area of the species

- env_occ:

  Data frame with two columns and 545 rows each representing
  environmental information bio1 and bio12 associated with the
  occurrence the species

- coords_m:

  Data frame with geographical coordinates of the accesitibility area
  (M).

- coords_occ:

  Data frame with geographical coordinates of the occurrence points.

## Examples

``` r
head(example_vicugna$env_occ)
#> # A tibble: 6 × 2
#>    bio1 bio12
#>   <dbl> <dbl>
#> 1  4.46   112
#> 2  5.74    98
#> 3  4.70    99
#> 4 14.7    127
#> 5 13.8    121
#> 6  1.59   109
```
