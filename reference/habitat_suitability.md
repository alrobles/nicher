# Tiled habitat-suitability map from an environmental terra stack

Evaluates the standardized multivariate-normal density of Jimenez et al.
(2022, Eq. 2) at every cell of an environmental
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html),
returning a one-layer raster of suitability values in \\(0, 1\]\\ (or,
optionally, log-suitability in \\(-\infty, 0\]\\).

## Usage

``` r
habitat_suitability(
  param,
  env,
  output = "",
  overwrite = FALSE,
  return_log = FALSE,
  threads = RcppParallel::defaultNumThreads(),
  wopt = list()
)
```

## Arguments

- param:

  Named list with components

  `mu`

  :   Numeric vector of length `p` — niche centroid.

  `Sigma`

  :   Symmetric positive-definite `p x p` matrix.

  The Cholesky factor of `Sigma` is computed once before the block loop.
  If [`chol()`](https://rdrr.io/r/base/chol.html) fails (rank-deficient
  `Sigma`), the function emits a
  [`warning`](https://rdrr.io/r/base/warning.html) and returns a raster
  of `NA`.

- env:

  A multi-layer
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html):
  one layer per environmental variable, in the same order as `param$mu`.
  Layer names are preserved on the input but are not used by this
  function — variable matching is positional. Use
  [`predict.nicher`](https://alrobles.github.io/nicher/reference/predict.nicher.md)
  when you need name-based matching.

- output:

  Character. File path for the output GeoTIFF. The empty string `""`
  (default) returns an in-memory `SpatRaster`.

- overwrite:

  Logical. If `TRUE` and `output` exists, it is overwritten. Default
  `FALSE`.

- return_log:

  Logical. If `FALSE` (default) the output is suitability \\S(x) \in (0,
  1\]\\. If `TRUE` the output is \\\log S(x) \le 0\\.

- threads:

  Integer. Number of parallel threads handed to the inner C++ kernel.
  Default
  [`RcppParallel::defaultNumThreads()`](https://rdrr.io/pkg/RcppParallel/man/setThreadOptions.html).

- wopt:

  List. Additional write options forwarded to
  [`writeStart`](https://rspatial.github.io/terra/reference/readwrite.html).

## Value

A one-layer `SpatRaster` named `"suitability"` (or `"log_suitability"`
when `return_log = TRUE`). Returned invisibly when `output` is
non-empty.

## Details

Memory is bounded by the largest block selected by terra's memory
manager: continental-scale rasters are processed without ever
materialising the full grid in R. The per-cell computation is run by
[`niche_suitability_cpp`](https://alrobles.github.io/nicher/reference/niche_suitability_cpp.md),
an RcppParallel kernel that takes raw pointers into terra's column-major
buffer.

The function follows the streaming I/O pattern from *xsdm-devel*'s
recipe 03:

1.  [`terra::readStart`](https://rspatial.github.io/terra/reference/readwrite.html)
    on `env`, paired with an `on.exit(terra::readStop)`;

2.  [`terra::writeStart`](https://rspatial.github.io/terra/reference/readwrite.html)
    on the output, paired with an `on.exit(terra::writeStop)`;

3.  for each block `i`: read each variable's tile via
    `terra::readValues(..., mat = TRUE)`, pack into a flat column-major
    numeric vector `(n_tile x p)`, mask `NA` cells, compact, call
    `niche_suitability_cpp`, scatter the result back, and
    [`terra::writeValues`](https://rspatial.github.io/terra/reference/readwrite.html).

## See also

[`niche_suitability_cpp`](https://alrobles.github.io/nicher/reference/niche_suitability_cpp.md),
[`predict.nicher`](https://alrobles.github.io/nicher/reference/predict.nicher.md).

## Examples

``` r
# \donttest{
  ## Hand-built (mu, Sigma):
  r1 <- terra::rast(nrows = 10, ncols = 10)
  r2 <- terra::rast(nrows = 10, ncols = 10)
  terra::values(r1) <- rnorm(100)
  terra::values(r2) <- rnorm(100)
  ex <- c(r1, r2)
  names(ex) <- c("bio1", "bio12")
  habitat_suitability(
    param = list(mu = c(0, 0), Sigma = diag(2)),
    env   = ex
  )
#> class       : SpatRaster
#> size        : 10, 10, 1  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84)
#> source(s)   : memory
#> name        : suitability
#> min value   :    0.008734
#> max value   :     0.97601
# }
```
