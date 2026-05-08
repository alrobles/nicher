# Multivariate Gaussian KDE with fixed Scott bandwidth

Computes a multivariate Gaussian kernel density estimate at the rows of
\`x\`, using \`data\` as the reference sample. Bandwidths are set
internally via Scott's rule-of-thumb (diagonal). For 2D, a manually
optimized version is used; for higher dimensions, an Eigen-based
vectorized version is used.

## Usage

``` r
kde_gaussian(x, data)
```

## Arguments

- x:

  Numeric matrix \`n_eval x p\` (evaluation points).

- data:

  Numeric matrix \`n_data x p\` (reference sample for KDE).

## Value

Numeric vector of length \`n_eval\` with the density estimates.

## Examples

``` r
x    <- as.matrix(example_env_occ_2d)
data <- as.matrix(example_env_m_2d)
dens <- kde_gaussian(x, data)
head(dens)
#> [1] 1.761649e-05 4.822381e-05 2.004299e-05 1.185014e-05 1.116945e-05
#> [6] 3.863498e-05
```
