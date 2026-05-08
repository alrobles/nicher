# Negative log-likelihood (inverse-probability-weighted normal, integrated C++)

Low-level C++ bridge for the inverse-probability-weighted (IPW) normal
niche model. Computes KDE weights and Mahalanobis distances entirely in
C++ to minimise R overhead. This is the workhorse called by
[`loglik_niche_math_ip_weighted`](https://alrobles.github.io/nicher/reference/loglik_niche_math_ip_weighted.md).

## Usage

``` r
loglik_niche_math_ip_weighted_integrated(
  theta,
  env_occ,
  env_m,
  eta = 1,
  neg = TRUE,
  den_idx = NULL,
  kde_idx = NULL,
  precomp_w_den = NULL
)
```

## Arguments

- theta:

  Numeric vector of unconstrained parameters (math scale): \\\[\mu,
  \log\sigma, v\]\\ of length \\2p + p(p-1)/2\\.

- env_occ:

  Data frame or matrix (\\n \times p\\) of environmental values at
  presence points.

- env_m:

  Data frame or matrix (\\n_m \times p\\) of background environmental
  values from the accessible area M.

- eta:

  Numeric scalar, shape parameter for the LKJ C-vine prior on the
  correlation matrix (default 1 = uniform over correlations).

- neg:

  Logical. If `TRUE` (default), returns the negative log-likelihood
  (suitable for minimisation).

- den_idx:

  Integer vector of 1-based row indices for the denominator subsample.
  If `NULL` (default), uses all rows of `env_m`.

- kde_idx:

  Integer vector of 1-based row indices for the KDE reference subsample.
  If `NULL`, uses all rows of `env_m`.

- precomp_w_den:

  Numeric vector of precomputed KDE weights for the denominator points
  (must match `den_idx` in length). If provided, KDE for the denominator
  is not recomputed.

## Value

Scalar numeric: the (negative) log-likelihood.

## Details

The objective function is described in detail in
[`loglik_niche_math_ip_weighted`](https://alrobles.github.io/nicher/reference/loglik_niche_math_ip_weighted.md).

This function was formerly called
`loglik_niche_math_kde_bias_corrected_integrated`.

## See also

\[loglik_niche_math_ip_weighted()\] for the user-facing wrapper with
automatic subsampling.

## Examples

``` r
# \donttest{
theta <- start_theta(example_env_occ_2d)
ll <- loglik_niche_math_ip_weighted_integrated(
  theta   = theta,
  env_occ = example_env_occ_2d,
  env_m   = example_env_m_2d,
  eta     = 1,
  neg     = TRUE
)
print(ll)
#> [1] 633.8731
# }
```
