# Optimize niche model log-likelihood with multi-start Sobol design

Runs multi-start optimization over a Sobol low-discrepancy sequence (via
pomp) of starting points covering the parameter space implied by
`env_occ` and `breadth`.

## Usage

``` r
optimize_niche(
  env_occ,
  env_m,
  num_starts = 100L,
  breadth = 0.1,
  likelihood = c("weighted", "ip_weighted", "presence_only", "skew_normal",
    "skew_normal_weighted", "skew_t", "skew_t_weighted"),
  grad = c("auto", "analytic", "central", "forward"),
  m_subsample = NULL,
  m_kde_subsample = NULL,
  seed = NULL,
  warm_start = TRUE,
  prior_log_sigma_lambda = 1,
  prior_log_sigma_center = NULL,
  prior_mu_lambda = 1,
  prior_mu_center = NULL,
  prior_alpha_lambda = 0.1,
  control = list(),
  verbose = FALSE,
  ...
)
```

## Arguments

- env_occ:

  Data frame of environmental values at presence points.

- env_m:

  Data frame of background environmental values. Required for
  `likelihood` in `"ip_weighted"`, `"weighted"`,
  `"skew_normal_weighted"`, or `"skew_t_weighted"`; ignored for
  `"presence_only"`, `"skew_normal"`, and `"skew_t"`.

- num_starts:

  Integer. Number of Sobol starting points.

- breadth:

  Numeric in (0, 0.5). Controls the quantile range used to define
  starting bounds for `mu` parameters. Default `0.1`.

- likelihood:

  One of `"weighted"` (default; paper Eq. 5 + ridge), `"ip_weighted"`
  (legacy KDE-bias-corrected formula), `"presence_only"`,
  `"skew_normal"` (presence-only multivariate skew-normal),
  `"skew_normal_weighted"` (paper Eq. 5 with skew-normal density + ridge
  prior on log sigma), `"skew_t"` (presence-only multivariate
  non-central skew-t via 32-node Gauss-Laguerre quadrature), or
  `"skew_t_weighted"` (paper Eq. 5 with NCST density + ridge prior).

- grad:

  Gradient strategy: `"auto"` (default) selects `"analytic"` for the
  Gaussian weighted models and `"central"` otherwise. Force one of
  `c("analytic", "central", "forward")` to override.

- m_subsample, m_kde_subsample:

  Optional integer or fraction in (0, 1\]. Resolved to
  `min(nrow(env_m), 10000)` when `NULL` (default).

- seed:

  Optional integer to make subsampling deterministic.

- warm_start:

  Logical. When `likelihood` is `"ip_weighted"`, `"weighted"`, or
  `"skew_normal_weighted"`, run a quick presence-only fit first and
  prepend its `theta` (padded with \\\alpha = 0\\ for the skew variants)
  as one extra starting point for the weighted multi-start. Cheap
  insurance against bad multistart luck on rough or multimodal weighted
  likelihoods (e.g. when `env_m` contains regions far from the
  occurrence cloud); does not replace the Sobol starts. Ignored for
  `likelihood = "presence_only"` and `"skew_normal"`. Default `TRUE`.

- prior_log_sigma_lambda:

  Numeric scalar (`>= 0`), strength of the ridge penalty on \\\log
  \sigma\\ used by `likelihood = "weighted"`, `"skew_normal_weighted"`,
  and `"skew_t_weighted"`. The penalty is \\\lambda \sum_k (\log
  \sigma_k - \log \hat\sigma_k)^2\\. The centre \\\log \hat\sigma_k\\
  defaults to the presence-only fit's \\\log \sigma_k\\ (shrinks the
  weighted fit toward the PO fit); see `prior_log_sigma_center`. Default
  `1.0` (weak). Larger values keep \\\sigma\\ closer to the PO scale;
  `0` reduces the model to pure ML weighted normal (Eq. 5) and is not
  recommended (subject to Patil & Ord 1976 drift).

- prior_log_sigma_center:

  Optional numeric vector of length `ncol(env_occ)`. `NULL` (default)
  anchors the centre at the presence-only fit's \\\log \sigma\\ when
  `warm_start = TRUE` (the PO fit is already run to seed the weighted
  multi-start); falls back to `log(apply(env_occ, 2, sd))` when
  `warm_start = FALSE`.

- prior_mu_lambda:

  Numeric scalar (`>= 0`), strength of a unit-free ridge penalty on
  \\\mu\\ used by the weighted families (`"weighted"`,
  `"skew_normal_weighted"`, `"skew_t_weighted"`). The penalty is
  \\\lambda\_\mu \sum_k \left((\mu_k - \hat\mu_k) /
  \hat\sigma_k\right)^2\\ with anchors \\\hat\mu_k\\ and
  \\\hat\sigma_k\\ from the PO fit (see `prior_mu_center`). Because the
  penalty is divided by the PO \\\sigma_k\\, it is scale-free across
  heterogeneous variables (e.g. bio1 in \\^\circ C\\ vs bio12 in mm):
  `prior_mu_lambda = 1` allows \\\mu\\ to drift ~1 PO standard deviation
  before the penalty pushes back. Default `1.0`. Set to `0` for
  unpenalized maximum likelihood.

- prior_mu_center:

  Optional numeric vector of length `ncol(env_occ)`. `NULL` (default)
  anchors the centre at the presence-only fit's \\\mu\\ when
  `warm_start = TRUE`; if `warm_start = FALSE` and
  `prior_mu_lambda > 0`, it must be supplied explicitly.

- prior_alpha_lambda:

  Numeric scalar (`>= 0`), strength of a ridge penalty \\\lambda\_\alpha
  \sum_k \alpha_k^2\\ on the skew vector \\\alpha\\, shrinking toward
  the Gaussian sub-model. Applies to both presence-only and weighted
  skew families (`"skew_normal"`, `"skew_normal_weighted"`, `"skew_t"`,
  `"skew_t_weighted"`). Fixes the well-known unbounded-MLE pathology of
  the skew-normal direct parameterization (Azzalini 1985, Pewsey 2000).
  Default `0.1` (mild). Set to `0` for unpenalized maximum likelihood.

- control:

  Named list of control parameters for
  [`ucminfcpp::ucminf_xptr()`](https://alrobles.github.io/ucminfcpp/reference/ucminf_xptr.html).
  Recognized entries:

  `grad`

  :   "central" (default)

  `gradstep`

  :   c(1e-6, 1e-8)

  `grtol`

  :   1e-4

  `xtol`

  :   1e-8

  `stepmax`

  :   5

  `maxeval`

  :   2000

- verbose:

  Logical. If `TRUE`, print per-start progress.

- ...:

  Additional arguments forwarded to the objective function (e.g. `eta`).

## Value

An object of class `"nicher"` (see
[`new_nicher`](https://alrobles.github.io/nicher/reference/new_nicher.md)).

## Details

The supported likelihood models are:

- `"weighted"` (default): paper-faithful weighted-normal model.
  Implements Eq. 5 of Jiménez & Soberón (2022, Ecological Modelling
  438:109982) – pure ML estimation of a weighted normal density on M –
  plus a weakly-informative ridge prior on \\\log \sigma\\ that prevents
  the well-known Patil & Ord (1976) \\\sigma \to \infty\\ drift in the
  un-regularised MLE. The ridge strength is controlled by
  `prior_log_sigma_lambda` (default `1.0`, weak); set it to `0` for pure
  paper Eq. 5 (not recommended on multimodal M).

- `"ip_weighted"`: inverse-probability-weighted (IPW) normal model.
  Assumes the Jimenez et al. (2019) DGP in which presences are drawn
  from \\f\\ restricted to \\M\\. Because background grid cells sample
  environmental space with density \\g(\mathbf{x})\\, a naive
  denominator sum estimates \\\int f\\g\\ rather than \\\int f\\. The
  IPW correction uses \\w = 1/\hat{g}\\ (the inverse of the background
  KDE) as importance-sampling weights (Eq. 8), yielding a
  Horvitz–Thompson estimator of the Lebesgue integral and recovering the
  fundamental niche on uniform environmental space. Empirically stable,
  with built-in self-regularisation against \\\sigma \to \infty\\ drift
  (no ridge prior required). The `"weighted"` model is preferred when
  observation intensity correlates with environmental density (e.g.
  opportunistic records); `"ip_weighted"` is preferred when sampling
  effort is approximately uniform across \\M\\ and the goal is to
  estimate intrinsic tolerances. See
  [`loglik_niche_math_ip_weighted`](https://alrobles.github.io/nicher/reference/loglik_niche_math_ip_weighted.md)
  for a detailed comparison and known limitations. Formerly called
  `"kde_bias_corrected"`; the old name is accepted as a deprecated
  alias.

- `"presence_only"`: model using only presence points, no background
  correction. Unimodal likelihood, useful as a sanity check or to seed
  the weighted multistart (see `warm_start`).

- `"skew_normal"`: presence-only fit of a multivariate skew-normal niche
  (Azzalini & Capitanio 1999, J. R. Stat. Soc. Ser. B 61(3): 579-602):
  \$\$S(x) \propto \phi_p(x - \mu; \Sigma) \\ \Phi\left( \sum_k \alpha_k
  \\ (x_k - \mu_k) / \sigma_k \right).\$\$ Adds a length-`p` skewness
  vector \\\alpha\\ to the existing \\(\mu, \Sigma)\\ parameters.
  \\\alpha = 0\\ recovers the symmetric Gaussian niche exactly.
  Sobol-start machinery samples \\\alpha_k\\ in \\\[-3, 3\]\\ (Azzalini
  & Capitanio 1999, Sec. 5).

- `"skew_normal_weighted"`: paper Eq. 5 with the SN density in place of
  the Gaussian, plus the same ridge prior on \\\log \sigma\\ used by
  `"weighted"`. Defaults inherit from `"weighted"`.

- `"skew_t"`: presence-only fit of the multivariate non-central skew-t
  (NCST) density (Branco & Dey 2001, J. Multivariate Analysis
  79(1):99-113): \$\$T = \mu + X \sqrt{r/Y}, \quad X \sim SN_p(0,
  \Sigma, \alpha), Y \sim \chi^2_r.\$\$ Adds a single degrees-of-freedom
  parameter `log_r` on top of the skew-normal layout. Smaller \\r\\ =
  heavier tails; \\r \to \infty\\ recovers `"skew_normal"`. The marginal
  density has no closed form, so it is computed by 32-node
  Gauss-Laguerre quadrature on the \\\chi^2_r\\ mixing variable.
  Sobol-start machinery samples `log_r` in \\\[\log 2, \log 100\]\\.

- `"skew_t_weighted"`: paper Eq. 5 with the NCST density in place of the
  Gaussian, plus the same ridge prior on \\\log \sigma\\ used by
  `"weighted"`. Defaults inherit from `"weighted"`.

For multimodal or long-tailed `env_m`, consider z-scoring `env_occ` and
`env_m` so all variables have comparable spread (e.g. z-score or
quantile-rank); the ridge prior is then uniform across axes.

## Optimization backend

`optimize_niche()` optimizes via
[`ucminfcpp::ucminf_xptr()`](https://alrobles.github.io/ucminfcpp/reference/ucminf_xptr.html)
with a compiled C++ objective built by `create_niche_obj_ptr()`. The
full `(theta -> mu, log sigma, v)` unpacking, `cvine_cholesky`,
log-likelihood, and gradient are evaluated in pure C++ with no
R-callback overhead. For the weighted likelihood, gradient mode
`"analytic"` uses closed-form derivatives over `mu` and `log sigma` and
central finite differences over the C-vine partial-correlation block.

## KDE sampling (weighted model)

The KDE that weights the M-background depends only on the environmental
data (not on `theta`), so the weights are computed once before
optimization. By default, the KDE reference (`m_kde_subsample`) and the
denominator subset (`m_subsample`) are both capped at 10\\000 distinct
background combinations – the maximum the model can usefully exploit
even for very large rasters. A floor of `max(500, 50 * 2^p)` rows is
used as a heuristic minimum representative sample (Silverman 1986; Wand
& Jones 1995); below this floor a
[`warning()`](https://rdrr.io/r/base/warning.html) is emitted.

## References

Jimenez, L., Soberon, J., Christen, J. A., & Soto, D. (2019). On the
problem of modeling a fundamental niche from occurrence data.
*Ecological Modelling*, 397, 109823.

Jimenez, L., & Soberon, J. (2022). Weighted-normal model for the
fundamental niche. *Ecological Modelling*, 438, 109982.

Azzalini, A., & Capitanio, A. (1999). Statistical applications of the
multivariate skew normal distribution. *J. R. Stat. Soc. B*, 61(3),
579–602.

Branco, M. D., & Dey, D. K. (2001). A general class of multivariate
skew-elliptical distributions. *J. Multivariate Anal.*, 79(1), 99–113.

## See also

\[loglik_niche_math_cpp()\], \[loglik_niche_math_presence_only()\],
\[loglik_niche_math_ip_weighted()\].

## Examples

``` r
# \donttest{
result <- optimize_niche(
  env_occ    = example_env_occ_2d,
  env_m      = example_env_m_2d,
  num_starts = 5L,
  breadth    = 0.1,
  likelihood = "weighted"
)
print(result)
#> -- nicher optimization result --
#>   Likelihood : weighted 
#>   Starts     : 6 (converged: 6 )
#>   Best loglik: -654.8488 
#>   Convergence: 1 
#>   eta        : 1 
#>   Variables  : bio1WH, bio12WH 
assess(result)
#> $flag
#> [1] "accepted_global"
#> 
#> $recommendation
#> [1] "Global optimum likely found: multiple converged starts agree closely."
#> 
#> $gap
#> [1] 1.932676e-12
#> 
#> $rel_gap
#> [1] 2.951332e-15
#> 
#> $n_converged
#> [1] 6
#> 
#> $best_loglik
#> [1] -654.8488
#> 
# }
```
