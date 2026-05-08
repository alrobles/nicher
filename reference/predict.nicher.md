# Habitat-suitability raster from a fitted `nicher` object

`predict` method for objects of class `"nicher"` returned by
[`optimize_niche`](https://alrobles.github.io/nicher/reference/optimize_niche.md).
Reconstructs the niche centroid \\\mu\\ and covariance \\\Sigma\\ from
the optimizer's `best$theta` (mu, log_sigma, v) parameterization, then
evaluates the standardized Gaussian suitability map of Jimenez et al.
(2022, Eq. 2) over an environmental
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html).

## Usage

``` r
# S3 method for class 'nicher'
predict(
  object,
  env,
  ...,
  return_log = FALSE,
  threads = RcppParallel::defaultNumThreads(),
  output = "",
  overwrite = FALSE,
  wopt = list()
)
```

## Arguments

- object:

  A `"nicher"` object returned by
  [`optimize_niche`](https://alrobles.github.io/nicher/reference/optimize_niche.md).

- env:

  A multi-layer
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html),
  one layer per environmental variable. When `object$var_names` is
  non-`NULL` (the default for fits produced by
  [`optimize_niche()`](https://alrobles.github.io/nicher/reference/optimize_niche.md)
  on a column-named `env_occ`), the layers of `env` are reordered to
  match `var_names` and any name not present in `env` produces an error.
  When `object$var_names` is `NULL` the layers are used in their
  existing order. This lets you call `predict` on a future-climate
  `SpatRaster` that contains the same variables in any order, including
  extra layers.

- ...:

  Currently ignored.

- return_log:

  Logical. `FALSE` (default) returns suitability in \\(0, 1\]\\; `TRUE`
  returns log-suitability in \\(-\infty, 0\]\\.

- threads:

  Integer. Threads for the inner C++ kernel. Default
  [`RcppParallel::defaultNumThreads()`](https://rdrr.io/pkg/RcppParallel/man/setThreadOptions.html).

- output:

  Character. File path for the output GeoTIFF. The empty string `""`
  (default) returns an in-memory raster.

- overwrite:

  Logical. Forwarded to
  [`habitat_suitability`](https://alrobles.github.io/nicher/reference/habitat_suitability.md).

- wopt:

  List. Forwarded to
  [`habitat_suitability`](https://alrobles.github.io/nicher/reference/habitat_suitability.md).

## Value

A one-layer `SpatRaster` named `"suitability"` (or `"log_suitability"`).
Returned invisibly when `output` is non-empty.

## Details

The same Gaussian suitability surface is used for every fitted family:
`"weighted"`, `"ip_weighted"`, `"presence_only"`, and the skew-normal /
skew-t variants all project the fitted \\(\mu, \Sigma)\\ geometry. For
skew fits, the skewness and tail parameters affect fitting but are not
used by this Gaussian projection.

Internally:

1.  \\\mu = \theta\_{1:p}\\.

2.  \\\sigma = \exp(\theta\_{p+1:2p})\\.

3.  C-vine partial-correlation block \\v =
    \theta\_{2p+1:\\length(\theta)}\\ is mapped to a correlation
    Cholesky factor via
    [`cvine_cholesky`](https://alrobles.github.io/nicher/reference/cvine_cholesky.md).

4.  \\L = \mathrm{diag}(\sigma)\\ L\_{corr}\\, \\\Sigma = L L^\top\\.

## See also

[`habitat_suitability`](https://alrobles.github.io/nicher/reference/habitat_suitability.md),
[`niche_suitability_cpp`](https://alrobles.github.io/nicher/reference/niche_suitability_cpp.md),
[`cvine_cholesky`](https://alrobles.github.io/nicher/reference/cvine_cholesky.md).

## Examples

``` r
# \donttest{
  fit <- optimize_niche(
    env_occ    = example_env_occ_2d,
    env_m      = example_env_m_2d,
    num_starts = 5L,
    likelihood = "weighted"
  )
  r1 <- terra::rast(nrows = 10, ncols = 10)
  r2 <- terra::rast(nrows = 10, ncols = 10)
  terra::values(r1) <- rnorm(100)
  terra::values(r2) <- rnorm(100)
  env <- c(r1, r2)
  names(env) <- colnames(example_env_occ_2d)
  suit <- predict(fit, env)
# }
```
