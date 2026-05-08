# Create a nicher optimization result object

Constructs an S3 object of class `"nicher"` from the outputs of
[`optimize_niche`](https://alrobles.github.io/nicher/reference/optimize_niche.md).

## Usage

``` r
new_nicher(solutions, best, likelihood, n_starts, var_names = NULL, eta = 1)
```

## Arguments

- solutions:

  A data frame with columns `start_id`, `loglik`, `convergence`, and
  `full_par` (list column).

- best:

  A list with the best solution: `theta`, `loglik`, `convergence`.

- likelihood:

  Character. One of `"ip_weighted"` or `"presence_only"`.

- n_starts:

  Integer. Total number of starting points used.

- var_names:

  Optional character vector of length `p` naming the environmental
  variables in the order assumed by `best$theta`. Stored on the object
  so
  [`predict.nicher`](https://alrobles.github.io/nicher/reference/predict.nicher.md)
  can reorder a future
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  by layer name.

- eta:

  Numeric scalar. The Beta-shape parameter passed to
  [`cvine_cholesky`](https://alrobles.github.io/nicher/reference/cvine_cholesky.md)
  during fitting (default `1.0`). Stored on the object so
  [`predict.nicher`](https://alrobles.github.io/nicher/reference/predict.nicher.md)
  reconstructs the same correlation matrix from `best$theta` that the
  optimizer used. Defaults to `1.0` (the LKJ-uniform prior), matching
  [`optimize_niche`](https://alrobles.github.io/nicher/reference/optimize_niche.md).

## Value

An object of class `"nicher"`.
