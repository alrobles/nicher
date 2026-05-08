# Negative log-likelihood (presence-only, math scale)

Computes the negative log-likelihood of a multivariate normal niche
model using only presence records, with no background correction
(Jimenez et al. 2019, Eq. 2 without the M-restricted denominator).

## Usage

``` r
loglik_niche_math_presence_only(theta, env_occ, eta = 1, neg = TRUE, ...)
```

## Arguments

- theta:

  Numeric vector of unconstrained parameters (math scale): \\\[\mu,
  \log\sigma, v\]\\ of length \\2p + p(p-1)/2\\.

- env_occ:

  Data frame or matrix (\\n \times p\\) of environmental values at
  presence points.

- eta:

  Numeric scalar, shape parameter for the LKJ C-vine prior on the
  correlation matrix (default 1 = uniform over correlations).

- neg:

  Logical. If `TRUE` (default), returns the negative log-likelihood
  (suitable for minimisation).

- ...:

  Additional arguments (ignored; for compatibility).

## Value

Scalar numeric: the (negative) log-likelihood.

## Details

The minimised objective is: \$\$ -\log\mathcal{L} =
\frac{n}{2}\log\|\Sigma\| + \frac{1}{2}\sum\_{i=1}^{n} (\mathbf{x}\_i -
\mu)^\top \Sigma^{-1}(\mathbf{x}\_i - \mu) \$\$ The \\(2\pi)\\
normalisation constant is dropped (does not affect the optimum).

Internally, \\\Sigma = L L^\top\\ is reconstructed from `theta` via
[`cvine_cholesky`](https://alrobles.github.io/nicher/reference/cvine_cholesky.md),
and the Mahalanobis distances are computed via a triangular solve
\\L^{-1}(\mathbf{x}\_i - \mu)\\.

## References

Jimenez, L., Soberon, J., Christen, J. A., & Soto, D. (2019). On the
problem of modeling a fundamental niche from occurrence data.
*Ecological Modelling*, 397, 109823.

## See also

\[optimize_niche()\] for multi-start fitting,
\[loglik_niche_math_cpp()\] for the M-restricted model.

## Examples

``` r
# \donttest{
theta <- start_theta(example_env_occ_2d)
ll <- loglik_niche_math_presence_only(theta, example_env_occ_2d)
# }
```
