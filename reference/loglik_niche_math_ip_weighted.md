# Negative log-likelihood (inverse-probability-weighted normal, math scale)

Computes the negative log-likelihood of the inverse-probability-weighted
(IPW) normal niche model. This model uses \\w = 1/\hat{g}\\ (the inverse
of the background KDE) as importance-sampling weights, applying a
Horvitz–Thompson correction so that the denominator approximates the
Lebesgue integral \\\int f(\mathbf{x})\\d\mathbf{x}\\ rather than the
\\g\\-weighted integral \\\int f\\g\\d\mathbf{x}\\:

## Usage

``` r
loglik_niche_math_ip_weighted(
  theta,
  env_occ,
  env_m,
  eta = 1,
  neg = TRUE,
  m_subsample = NULL,
  m_kde_subsample = NULL,
  seed = NULL,
  den_idx = NULL,
  kde_idx = NULL,
  precomp_w_den = NULL,
  ...
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

  Logical. If `TRUE` (default), returns negative log-likelihood.

- m_subsample:

  Optional integer or fraction. If `den_idx` is not given, this defines
  the number (or fraction) of rows of `env_m` used for the denominator
  subsample.

- m_kde_subsample:

  Optional integer or fraction. If `kde_idx` is not given, this defines
  the KDE reference subsample size (or fraction).

- seed:

  Optional integer seed for reproducible subsampling.

- den_idx:

  Optional integer vector of 1-based row indices for the denominator
  subset. If provided, no new denominator indices are generated.

- kde_idx:

  Optional integer vector of 1-based row indices for the KDE reference
  subset. If provided, no new KDE indices are generated.

- precomp_w_den:

  Optional numeric vector of precomputed denominator KDE weights. Must
  match the size of `den_idx`. If provided, KDE for the denominator is
  not recomputed.

- ...:

  Additional arguments (ignored; provided for compatibility).

## Value

A scalar numeric value containing the (negative) log-likelihood.

## Details

\$\$ -\log\mathcal{L} = \frac{1}{2}\sum_i q_i + \sum_i
\log\hat{g}(\mathbf{x}\_i) + n\cdot\mathrm{logsumexp}\_j\\\left\[
-\tfrac{1}{2}q_j - \log\hat{g}(\mathbf{y}\_j) \right\] \$\$

where \\q_i\\ is the Mahalanobis distance at occurrence point
\\\mathbf{x}\_i\\, and the sum over \\j\\ runs over background points
\\\mathbf{y}\_j \in M\\.

## Data-generating process and motivation

The IPW model is designed for the data-generating process (DGP) of
Jimenez et al. (2019), which assumes that presence records are drawn
from the fundamental niche density \\f(\mathbf{x};\mu,\Sigma)\\
restricted to the environments that actually exist in the accessible
area \\M\\:

\$\$ p(\mathbf{x}\_i \mid \text{presence in } M) =
\frac{f(\mathbf{x}\_i)}{\int_M f(\mathbf{x})\\d\mathbf{x}}. \$\$

The computational challenge is approximating the denominator \\\int_M
f\\d\mathbf{x}\\. The background sample \\\\\mathbf{y}\_j\\\\ is drawn
from geographic grid cells whose density in environmental space is
\\g(\mathbf{x})\\. A naive sum \\\frac{1}{K}\sum_j f(\mathbf{y}\_j)\\
therefore estimates \\\int f\\g\\d\mathbf{x}\\, not \\\int
f\\d\mathbf{x}\\. Dividing each term by \\\hat{g}(\mathbf{y}\_j)\\
corrects the sampling bias via importance sampling:

\$\$ \int_M f(\mathbf{x})\\d\mathbf{x} = \int_M
\frac{f(\mathbf{x})}{g(\mathbf{x})}\\g(\mathbf{x})\\d\mathbf{x} \approx
\frac{1}{K}\sum\_{j=1}^{K}
\frac{f(\mathbf{y}\_j)}{\hat{g}(\mathbf{y}\_j)}. \$\$

This is the classical Horvitz–Thompson estimator (Horvitz & Thompson
1952) applied to niche modelling. The resulting likelihood estimates the
fundamental niche on uniform environmental space (Lebesgue measure),
removing the distortion that arises when common environments in \\M\\
dominate the denominator.

## Comparison with the `"weighted"` model

The `"weighted"` model (Jimenez & Soberon 2022, Eq. 5/8) uses \\w =
\hat{g}\\, which is the correct likelihood under an alternative DGP
where presence intensity is proportional to \\f \cdot g\\ (the
use-availability or resource-selection function model; Warton & Shepherd
2010). The two models answer different ecological questions:

|  |  |  |
|----|----|----|
| **Model** | **DGP assumption** | **Estimates** |
| `"weighted"` | \\\lambda \propto f \cdot g\\ | Niche under M-biased observation process |
| `"ip_weighted"` | \\f\\ restricted to \\M\\ | Fundamental niche on uniform E-space |

Neither model is universally superior. `"weighted"` is appropriate when
observation probability correlates with environmental density (e.g.
opportunistic citizen-science records). `"ip_weighted"` is appropriate
when sampling effort is approximately uniform across \\M\\ and the goal
is to recover the species' intrinsic tolerances independent of which
environments happen to be common.

## Known limitations

1.  **Importance-sampling variance.** When \\\hat{g}(\mathbf{y}\_j)
    \approx 0\\ (e.g. between disconnected environmental patches in a
    multimodal \\M\\), the weights \\1/\hat{g}\\ can be very large and a
    handful of background points may dominate the denominator. The
    effective sample size (ESS) should be monitored; low ESS indicates
    unreliable estimates.

2.  **KDE bandwidth dependence.** \\\hat{g}\\ is computed once using
    Scott's rule bandwidth. Over- or under-smoothing can distort the
    weights. An adaptive or cross-validated bandwidth would improve
    robustness.

3.  **Self-regularisation.** The IPW denominator diverges when \\\sigma
    \to \infty\\ because rare-environment cells contribute \\1/\hat{g}
    \gg 1\\. This acts as a built-in penalty against unrealistically
    broad niches (no explicit ridge prior is needed), but the effective
    constraint depends on the tail behaviour of \\\hat{g}\\ and can be
    noisy.

Accepts explicit subsampling indices and precomputed KDE weights for
high-performance workflows. Internally delegates to
[`loglik_niche_math_ip_weighted_integrated`](https://alrobles.github.io/nicher/reference/loglik_niche_math_ip_weighted_integrated.md).

This function was formerly called
`loglik_niche_math_kde_bias_corrected`.

## References

Jimenez, L., Soberon, J., Christen, J. A., & Soto, D. (2019). On the
problem of modeling a fundamental niche from occurrence data.
*Ecological Modelling*, 397, 109823.

Jimenez, L., & Soberon, J. (2022). Weighted-normal model for the
fundamental niche. *Ecological Modelling*, 438, 109982.

Horvitz, D. G., & Thompson, D. J. (1952). A generalization of sampling
without replacement from a finite universe. *J. Amer. Statist. Assoc.*,
47(260), 663–685.

Warton, D. I., & Shepherd, L. C. (2010). Poisson point process models
solve the “pseudo-absence problem” for presence-only data in ecology.
*Ann. Appl. Stat.*, 4(3), 1383–1402.

## See also

\[optimize_niche()\] for multi-start fitting,
\[loglik_niche_math_cpp()\] for the M-restricted Gaussian model.

## Examples

``` r
# \donttest{
den_idx <- sample.int(nrow(example_env_m_2d), 2000)
kde_idx <- sample.int(nrow(example_env_m_2d), 5000)
pre_w <- kde_gaussian(
  example_env_m_2d[den_idx, ],
  example_env_m_2d[kde_idx, ]
)

loglik_niche_math_ip_weighted(
  theta   = start_theta(example_env_occ_2d),
  env_occ = example_env_occ_2d,
  env_m = example_env_m_2d,
  den_idx = den_idx,
  kde_idx = kde_idx,
  precomp_w_den = pre_w
)
#> [1] 521.1715
# }
```
