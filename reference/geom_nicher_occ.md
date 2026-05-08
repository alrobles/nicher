# Occurrence-points layer

A self-contained `ggplot2` layer that draws an occurrence cloud
`env_occ` in 2-D environmental space.

## Usage

``` r
geom_nicher_occ(
  env_occ,
  var_names = NULL,
  size = 1.2,
  alpha = 1,
  colour = "black",
  shape = 16,
  ...
)
```

## Arguments

- env_occ:

  Matrix or data frame with at least two columns. Only the first two
  columns are used.

- var_names:

  Optional character vector of length 2 used to name the columns of the
  layer's internal data frame.

- size, alpha, colour, shape:

  Aesthetic parameters.

- ...:

  Additional fixed parameters passed to `geom_point`.

## Value

A `ggplot2` layer.

## Examples

``` r
# \donttest{
  library(ggplot2)
  ggplot() +
    geom_nicher_background(example_env_m_2d) +
    geom_nicher_occ(example_env_occ_2d)

# }
```
