# Negative log likelihood of an ellipsoid corrected with environmental combinations which come from the area of study (M)

Negative log likelihood of an ellipsoid corrected with environmental
combinations which come from the area of study (M)

## Usage

``` r
loglik_niche(mu, s_mat, env_occ, env_m, neg = TRUE)
```

## Arguments

- mu:

  A vector mu of parameters

- s_mat:

  The covariance matrix from environmental data frame

- env_occ:

  A data.frame containing the original sample of environmental
  combinations that correspond to presences

- env_m:

  A data.frame containing a second random sample of environmental

- neg:

  Logical. Default TRUE, returns the negative of the likelihood

## Value

A negative log likelihood value

## Examples

``` r
loglik_niche(
  mu = example_mu_vec,
  s_mat = example_s_mat,
  env_occ = example_env_occ_2d,
  env_m = example_env_m_2d
)
#> [1] 661.8376
# Example with log of the likelihood
loglik_niche(
  mu = example_mu_vec,
  s_mat = example_s_mat,
  env_occ = example_env_occ_2d,
  env_m = example_env_m_2d,
  neg = FALSE
)
#> [1] -661.8376
```
