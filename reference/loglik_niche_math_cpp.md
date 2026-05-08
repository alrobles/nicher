# Negative log-likelihood (M-restricted, math scale, Cholesky version)

Computes the negative log-likelihood of the multivariate normal niche
model restricted to the set of existing environments \\\mathbf{E}(t;G)\\
(Jimenez et al. 2019, Eq. 3–5). The density at an occurrence point is
normalised by the sum of the density over all background points in M:

## Usage

``` r
loglik_niche_math_cpp(theta, env_occ, env_m, eta = 1, neg = TRUE, ...)
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

- ...:

  Additional arguments (ignored; for compatibility).

## Value

A scalar numeric value: the (negative) log-likelihood.

## Details

\$\$ -\log\mathcal{L} = \frac{1}{2}\sum\_{i=1}^{n} (\mathbf{x}\_i -
\mu)^\top \Sigma^{-1}(\mathbf{x}\_i - \mu) + n \cdot \log\\\left(
\sum\_{\mathbf{y} \in M} \exp\\\left\[-\tfrac{1}{2} (\mathbf{y} -
\mu)^\top \Sigma^{-1}(\mathbf{y} - \mu) \right\] \right) \$\$

The \\\|\Sigma\|^{-1/2}\\ factor cancels between numerator and
denominator (Eq. 3), so it does not appear in the objective. The
logsumexp is computed with the max-shift trick for numerical stability.

## References

Jimenez, L., Soberon, J., Christen, J. A., & Soto, D. (2019). On the
problem of modeling a fundamental niche from occurrence data.
*Ecological Modelling*, 397, 109823.

## See also

\[optimize_niche()\] for multi-start fitting,
\[loglik_niche_math_presence_only()\] for the unconstrained (no M)
version.

## Examples

``` r
# \donttest{
theta <- start_theta(example_env_occ_2d)
ll <- loglik_niche_math_cpp(theta,
  env_occ = example_env_occ_2d,
  env_m = example_env_m_2d,
  eta = 1, neg = TRUE
)
print(ll)
#> [1] 662.9334
# }
```
