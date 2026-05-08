# Habitat-suitability kernel (parallel, zero-copy)

Evaluates the standardized multivariate-normal density \\S(x) =
\exp(-\tfrac12 (x-\mu)^\top \Sigma^{-1} (x-\mu))\\ for each row of a
column-major flat environmental buffer. Uses RcppParallel's
`parallelFor` over locations.

## Usage

``` r
niche_suitability_cpp(
  env_dat_vec,
  env_dat_dims,
  mu,
  L_inv,
  return_log = FALSE,
  num_threads = 0L
)
```

## Arguments

- env_dat_vec:

  Numeric vector. Flat column-major buffer of length `n_loc * p`; entry
  for location `l` variable `k` sits at index `l + n_loc * k`.

- env_dat_dims:

  Integer vector `c(n_loc, p)`.

- mu:

  Numeric vector of length `p`: niche centroid.

- L_inv:

  Numeric matrix `p x p`: inverse of the lower Cholesky factor of
  \\\Sigma\\. Precomputed R-side once.

- return_log:

  Logical. If `FALSE` (default) returns suitability in \\(0, 1\]\\; if
  `TRUE` returns \\\log S(x) \le 0\\.

- num_threads:

  Integer. `0` (default) leaves `RcppParallel`'s global thread state
  unchanged.

## Value

Numeric vector of length `n_loc`.
