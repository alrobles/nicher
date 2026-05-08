# Fit presence-only Gaussian niche model

Fits a multivariate normal niche model using only presence records, with
no background correction (Jimenez et al. 2019, Eq. 2). Uses the compiled
C++ backend via an external pointer for fast evaluation.

## Usage

``` r
niche_presence_only(occ, eta = 1, start = NULL, ...)
```

## Arguments

- occ:

  Numeric matrix (\\n \times p\\) of environmental values at presence
  points.

- eta:

  Numeric scalar, shape parameter for the LKJ C-vine prior on the
  correlation matrix (default 1 = uniform over correlations).

- start:

  A numeric vector (single-start) or list of numeric vectors
  (multi-start). Use \[start_theta()\] or \[start_theta_multiple()\] to
  generate starting values.

- ...:

  Ignored. Present for wrapper compatibility.

## Value

A list with components:

- `theta`:

  Best parameter vector.

- `value`:

  Negative log-likelihood at the optimum.

- `conv`:

  Convergence code (0 = success).

- `all_results`:

  Data frame of all starts (multi-start only).

## Details

The minimised objective is the negative log-likelihood of the
unconstrained multivariate normal: \$\$ -\log\mathcal{L} =
\frac{n}{2}\log\|\Sigma\| + \frac{1}{2}\sum\_{i=1}^{n} q_i \$\$ where
\\q_i = (\mathbf{x}\_i - \mu)^\top \Sigma^{-1} (\mathbf{x}\_i - \mu)\\.

Supports both single-start (`start` is a numeric vector) and multi-start
(`start` is a list of numeric vectors) optimisation.

## References

Jimenez, L., Soberon, J., Christen, J. A., & Soto, D. (2019). On the
problem of modeling a fundamental niche from occurrence data.
*Ecological Modelling*, 397, 109823.

## See also

\[optimize_niche()\] for the unified fitting interface,
\[loglik_niche_math_presence_only()\] for the R-level log-likelihood.

## Examples

``` r
occ    <- as.matrix(example_env_occ_3d)
theta0 <- start_theta(example_env_occ_3d)
res    <- niche_presence_only(occ = occ, start = theta0)
res$value
#> [1] 988.4879
```
