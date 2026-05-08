# Fit inverse-probability-weighted (IPW) normal niche model

Fits the inverse-probability-weighted (IPW) normal niche model using a
compiled C++ backend via an external pointer for fast evaluation.

## Usage

``` r
niche_ip_weighted(
  occ,
  M,
  den_idx,
  kde_idx,
  precomp_w_den,
  eta = 1,
  start = NULL,
  ...
)
```

## Arguments

- occ:

  Numeric matrix (\\n \times p\\) of environmental values at presence
  points.

- M:

  Numeric matrix (\\n_m \times p\\) of environmental values from the
  accessible area \\M\\. Must have the same columns as `occ`.

- den_idx:

  Integer vector of 1-based row indices selecting the denominator
  subset. Must have the same length as `precomp_w_den`.

- kde_idx:

  Integer vector of 1-based row indices selecting the KDE reference
  subset.

- precomp_w_den:

  Numeric vector of precomputed KDE weights matching `den_idx` in
  length.

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

This model assumes the data-generating process of Jimenez et al. (2019):
presence records are drawn from the fundamental niche density
\\f(\mathbf{x};\mu,\Sigma)\\ restricted to the environments available in
\\M\\. Because background grid cells are not uniformly distributed in
environmental space (their density is \\g(\mathbf{x})\\), a naive sum
over the background would estimate \\\int f\\g\\d\mathbf{x}\\ rather
than the required \\\int f\\d\mathbf{x}\\. The IPW correction divides
each background term by \\\hat{g}\\ (a kernel density estimate of the
environmental density of \\M\\), yielding a Horvitz–Thompson estimator
(Horvitz & Thompson 1952) of the Lebesgue integral.

The objective is the negative log-likelihood: \$\$ -\log\mathcal{L} =
\frac{1}{2}\sum_i q_i + \sum_i \log\hat{g}(\mathbf{x}\_i) +
n\cdot\mathrm{logsumexp}\_j\\\left\[ -\tfrac{1}{2}q_j -
\log\hat{g}(\mathbf{y}\_j) \right\] \$\$

The `"weighted"` model (Jimenez & Soberon 2022, Eq. 5/8) uses \\w =
\hat{g}\\ instead, which is the correct likelihood under the alternative
use-availability DGP where presence intensity is proportional to \\f
\cdot g\\. See
[`loglik_niche_math_ip_weighted`](https://alrobles.github.io/nicher/reference/loglik_niche_math_ip_weighted.md)
for a detailed comparison of both DGP assumptions and known limitations.

This function was formerly called `niche_kde_bias_corrected`.

KDE weights must be precomputed by the user and passed via
`precomp_w_den`. No KDE is recomputed inside the optimiser.

Supports both single-start (`start` is a numeric vector) and multi-start
(`start` is a list of numeric vectors) optimisation.

## References

Jimenez, L., Soberon, J., Christen, J. A., & Soto, D. (2019). On the
problem of modeling a fundamental niche from occurrence data.
*Ecological Modelling*, 397, 109823.

Jimenez, L., & Soberon, J. (2022). Weighted-normal model for the
fundamental niche. *Ecological Modelling*, 438, 109982.

Horvitz, D. G., & Thompson, D. J. (1952). A generalization of sampling
without replacement from a finite universe. *J. Amer. Statist. Assoc.*,
47(260), 663–685.

## See also

\[optimize_niche()\] for the unified fitting interface,
\[loglik_niche_math_ip_weighted()\] for the R-level log-likelihood
(includes full DGP derivation and known limitations).

## Examples

``` r
# \donttest{
occ <- as.matrix(example_env_occ_3d)
M   <- as.matrix(example_env_m_3d)
set.seed(1)
den_idx <- sample.int(nrow(M), 300L)
kde_idx <- sample.int(nrow(M), 600L)
w_den   <- kde_gaussian(M[den_idx, ], M[kde_idx, ])
theta0  <- start_theta(example_env_occ_3d)
res <- niche_ip_weighted(occ, M, den_idx, kde_idx, w_den, start = theta0)
res$value
#> [1] 347.8219
# }
```
