# Background environment layer

A self-contained `ggplot2` layer that draws the environmental background
`env_m` as a faint point cloud in 2-D environmental space. Designed to
be composed with
[`geom_nicher_occ`](https://alrobles.github.io/nicher/reference/geom_nicher_occ.md)
and
[`geom_nicher_ellipse`](https://alrobles.github.io/nicher/reference/geom_nicher_ellipse.md)
via the `+` operator.

## Usage

``` r
geom_nicher_background(
  env_m,
  var_names = NULL,
  size = 0.3,
  alpha = 0.25,
  colour = "grey60",
  ...
)
```

## Arguments

- env_m:

  Matrix or data frame with at least two columns. Only the first two
  columns are used (positional indexing).

- var_names:

  Optional character vector of length 2 used to name the columns of the
  layer's internal data frame. Defaults to `c("x1", "x2")`.

- size, alpha, colour:

  Aesthetic parameters; defaults are tuned for a discreet background
  layer.

- ...:

  Additional fixed parameters passed to the underlying
  [`ggplot2::geom_point`](https://ggplot2.tidyverse.org/reference/geom_point.html).

## Value

A `ggplot2` layer.

## Examples

``` r
# \donttest{
  library(ggplot2)
  ggplot() + geom_nicher_background(example_env_m_2d)

# }
```
