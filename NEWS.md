# nicher 3.2.0

## Information criteria + side-by-side model comparison (new feature, non-breaking)

The 7-likelihood + multi-knob landscape introduced in 3.0.x / 3.1.0
made it easy to fit many candidate models. This release adds the
standard frequentist tooling to compare them.

### New S3 methods

- `logLik()` -- returns the **un-penalised** log-likelihood evaluated
  at the converged `theta`, with `df = length(theta)` and `nobs =
  nrow(env_occ)`. Important: this is **not** the same as
  `fit$best$loglik` for penalised fits (`fit$best$loglik` is the
  optimiser's objective, which embeds the ridge penalty). AIC and BIC
  are defined on the bare data likelihood, so `logLik(fit)` is what
  `AIC()` / `BIC()` consume.
- `nobs()` -- the sample size used by BIC (number of occurrence rows).
- `AIC()`, `BIC()` -- standard formulas (`-2 ll + 2 k` and
  `-2 ll + k log n`). Pass multiple fits to get a one-row-per-fit
  data.frame.

### New `compare_nicher()` helper

```
compare_nicher(weighted = fit_w,
               skew_normal = fit_skn,
               skew_t      = fit_skt)
```

Returns a `data.frame(model, likelihood, loglik, df, nobs, AIC, dAIC,
BIC, dBIC, weight_AIC, convergence)` ranked best-first by AIC (or BIC,
or raw `loglik`). The Akaike weights `weight_AIC` quantify each
model's posterior probability under a flat model prior:
`exp(-d_i / 2) / sum_j exp(-d_j / 2)`.

`compare_nicher()` validates that all fits used the **same input
data**: identical `env_occ` fingerprint and (when applicable)
identical `env_m` fingerprint. Mismatches raise a hard error rather
than silently producing meaningless IC values. Mixing fits that use
`env_m` (`weighted`, `skew_normal_weighted`, `skew_t_weighted`,
`kde_bias_corrected`) with fits that do not (`presence_only`,
`skew_normal`, `skew_t`) is allowed but emits a warning -- weighted
likelihoods are normalised by an `env_m` integral and so live on a
different scale than presence-only likelihoods.

### Effective degrees of freedom

`df` is reported naively as `length(theta)`; this ignores the
shrinkage induced by the ridge penalties added in 3.1.0
(`prior_mu_lambda`, `prior_log_sigma_lambda`, `prior_alpha_lambda`).
For weak penalties (`lambda <= 1`) the bias is small. For stronger
penalties the effective df is lower than `length(theta)`, so AIC/BIC
will under-favour the more-regularised model. Quantifying effective
df via the trace of the influence matrix is left to a future release.

### Storage on the nicher object

Each `optimize_niche()` call now also stores on `fit$best`:

- `loglik_unpenalised` -- bare data log-likelihood at the optimum.
- `nobs` -- sample size.
- `env_occ_fingerprint`, `env_m_fingerprint` -- cheap deterministic
  hashes (nrow, ncol, colnames, column means, column SDs) used by
  `compare_nicher()` to detect non-comparable inputs.

These are additive; older `nicher` objects (`< 3.2.0`) cannot be
passed to `logLik()` / `AIC()` / `BIC()` and will raise a clear
"refit with current version" error.

# nicher 3.1.0

## PO-anchored penalised MLE (new feature, non-breaking)

The weighted families now shrink their fits toward the presence-only (PO)
baseline via three optional ridge penalties layered on the negative
log-likelihood. This operationalises the package's philosophy that the
accessible-environment term `env_m` should *correct* the PO niche, not
replace it: without a penalty, the weighted MLE for `mu` can drift many
standard deviations away from the PO centre on real data (on the
bundled `example_vicugna` dataset, unregularised `weighted` puts `mu2`
at 680 vs. the PO value of 215).

This is **penalised maximum likelihood** (Tikhonov / ridge), not a
Bayesian posterior: the optimiser still minimises `-log_lik +
penalty`; the penalty simply prefers fits that stay close to the PO
anchor.

### New `optimize_niche()` arguments

- `prior_mu_lambda`  (default `1.0`) -- unit-free ridge
  `sum_k ((mu_k - mu_PO_k) / sigma_PO_k)^2`. Because the penalty is
  scaled by the PO standard deviation, `lambda = 1` means "allow `mu`
  to drift ~1 PO sigma before the penalty pushes back". Scale-free
  across heterogeneous variables (e.g. `bio1` in degrees C vs. `bio12`
  in mm).
- `prior_alpha_lambda`  (default `0.1`) -- ridge `sum_k alpha_k^2`
  on the skew vector, shrinking toward the Gaussian sub-model. Fixes
  the well-known unbounded-MLE pathology of the skew-normal direct
  parameterization (Azzalini 1985, Pewsey 2000). Applies to **all**
  skew families, both presence-only (`skew_normal`, `skew_t`) and
  weighted (`skew_normal_weighted`, `skew_t_weighted`). Without this
  penalty, `example_vicugna` produced `alpha = (66 229, 1 287 048)`,
  clearly a degenerate optimum; with `lambda = 0.1` the fit lands at
  `alpha = (0.66, 13)`.
- `prior_mu_center`  (default `NULL`) -- centre of the `mu` ridge.
  When `NULL`, it is anchored at the PO fit's `mu` (reuses the
  warm-start optimisation that was already run to seed the weighted
  multi-start).

### Changed defaults

- `prior_log_sigma_lambda` and `prior_alpha_lambda` now **default to
  nonzero values** (`1.0` and `0.1` respectively) on their applicable
  families. To recover the pure-MLE behaviour of nicher 3.0.x set
  all three `prior_*_lambda = 0`.
- `prior_log_sigma_center` now defaults to the PO fit's `log(sigma)`
  (via warm-start) instead of `log(apply(env_occ, 2, sd))`. This is
  tighter and more honest -- the two agree up to numerical precision
  for the default `breadth`, but can differ on short or multimodal
  occurrence clouds.
- `warm_start = TRUE` (already the default) is now required for any
  weighted-family call that relies on the default anchors. Call with
  `warm_start = FALSE` and you must either supply `prior_mu_center`
  explicitly or set `prior_mu_lambda = 0`.

### Migration

- Fits from nicher 3.0.x with `prior_log_sigma_lambda = 1` (the 3.0
  default) and default `prior_log_sigma_center = log(sd(env_occ))`
  **will change** at nicher 3.1.0 because the centre now tracks the
  PO fit. To reproduce the 3.0 numerics exactly, pass
  `prior_log_sigma_center = log(apply(env_occ, 2, sd))`,
  `prior_mu_lambda = 0`, `prior_alpha_lambda = 0`.

### Bibliography

- Azzalini, A. (1985). A class of distributions which includes the
  normal ones. *Scand. J. Statist.* 12, 171-178.
- Pewsey, A. (2000). Problems of inference for Azzalini's skewnormal
  distribution. *J. Applied Statistics* 27(7), 859-870.

# nicher 3.0.0

## Breaking changes

### Hard rename of the weighted likelihood family

The `likelihood` strings accepted by `optimize_niche()` have been
**renamed and reordered** so that the **paper-faithful formula
(Jiménez & Soberón 2022, Eq. 5) is now the default**. This is a
breaking API change with no soft-deprecation alias: code passing
`likelihood = "weighted"` to nicher 2.x will silently get a different
formula in nicher 3.0, and code passing
`likelihood = "weighted_penalized"` will error.

| nicher 2.x                          | nicher 3.0                          | What it computes                                                          |
|-------------------------------------|-------------------------------------|---------------------------------------------------------------------------|
| `likelihood = "weighted"` (default) | `likelihood = "kde_bias_corrected"` | Legacy KDE-bias-corrected formula (production formula in 2.x)             |
| `likelihood = "weighted_penalized"` | `likelihood = "weighted"` (default) | Paper-faithful weighted normal (Eq. 5) + ridge prior on `log sigma`       |

To preserve the previous behaviour exactly, change
`likelihood = "weighted"` to `likelihood = "kde_bias_corrected"`.
To opt into the paper-faithful default, change
`likelihood = "weighted_penalized"` to `likelihood = "weighted"`.

### Why the rename?

The 2.x `likelihood = "weighted"` did *not* implement Eq. 5 of
Jiménez & Soberón (2022) -- it computed a sampling-bias-corrected
presence-only likelihood with a sign-flipped denominator term
(`exp(-q2/2) / w(y_j)` instead of `w(y_j) * exp(-q2/2)`). Both
formulas have valid statistical interpretations (the 2.x formula is
empirically more stable; the paper formula is the published one),
but it is misleading to call the non-paper formula `"weighted"`.
The paper formula takes the canonical name in 3.0; the legacy
formula gets a descriptive name reflecting what it actually
computes.

## New features

### Paper-faithful weighted-likelihood (`likelihood = "weighted"`)

* `likelihood = "weighted"` is now the default for `optimize_niche()`.
  Implements Eq. 5 of Jiménez & Soberón (2022, *Ecological Modelling*
  **438**: 109982) exactly -- pure ML estimation of a weighted normal
  density on M -- plus a weakly-informative ridge penalty on
  `log sigma`:
  ```
  penalty = lambda * sum_k (log sigma_k - log sigma_hat_k)^2
  ```
  with `log sigma_hat_k = log(sd(env_occ[, k]))`. The penalty
  prevents the well-known Patil & Ord (1976) `sigma -> infinity` drift
  that affects the un-regularised paper MLE on datasets where the
  background `env_m` is multimodal or has long tails relative to the
  occurrence cloud (e.g. `bio12` in the bundled `example_vicugna`
  dataset). Set `prior_log_sigma_lambda = 0` to disable the ridge and
  recover pure paper Eq. 5 (not recommended on multimodal M).
* The legacy KDE-bias-corrected formula remains available under
  `likelihood = "kde_bias_corrected"`. It is empirically stable and is
  what nicher 2.x users were running by default.
* Two new arguments on `optimize_niche()`:
    * `prior_log_sigma_lambda` -- ridge strength (default `1.0`,
      weakly informative). `0` disables the penalty (not recommended).
    * `prior_log_sigma_center` -- prior centre on `log sigma`; defaults
      to `log(apply(env_occ, 2, sd))`.
* For best results with `"weighted"` on real-world bioclim
  data, also rescale `env_occ` / `env_m` so all variables have
  comparable spread (e.g. z-score or quantile-rank both before
  fitting). The penalty is in `log sigma` units and a rescaling
  brings the per-axis scales into the same range.

### Warm-start for the weighted family

* New `warm_start = TRUE` argument on `optimize_niche()` (default
  enabled). When `likelihood` is `"kde_bias_corrected"`, `"weighted"`,
  or `"skew_normal_weighted"`, a quick presence-only fit is run first
  (~20 starts) and its `theta` is prepended as one extra starting
  point for the weighted multi-start. The presence-only likelihood is
  unimodal and converges cleanly, so this is cheap insurance against
  bad multistart luck on rough or multimodal weighted likelihoods --
  it does **not** replace the Sobol starts. For `skew_normal_weighted`
  the warm-start `theta` is padded with `alpha = 0`, which is the
  Gaussian (symmetric) point in the skew-normal parameter space.
* Failure to converge in the presence-only sub-fit is downgraded to
  a `warning()`; the weighted multi-start continues with the Sobol
  starts only.

### Skew-t niche likelihoods (`skew_t`, `skew_t_weighted`)

* Two new heavy-tailed likelihood families on `optimize_niche()` that
  fit a multivariate non-central skew-t niche (Branco & Dey 2001,
  *J. Multivariate Anal.* **79**(1), 99--113):
  ```
  T = mu + X * sqrt(r / Y),   X ~ SN_p(0, Sigma, alpha),  Y ~ chi^2_r
  ```
  Adds a single positive degrees-of-freedom parameter `r > 0` to the
  skew-normal parameter set (parameterised on `log r` for unconstrained
  optimisation). As `r -> Inf` the family reduces to the skew-normal
  exactly. The skew-t niche is useful when the environmental
  conditions at occurrences have heavier tails than the Gaussian or
  skew-normal allow -- a common situation when occurrences span large
  geographic gradients.
* The NCST density has no closed form; the per-point integral over the
  chi-squared mixing variable is approximated by 32-node standard
  Gauss--Laguerre quadrature (`Q = 32`, accuracy ~1e-8 over `r` in
  `[2, 100]`). Cost per likelihood evaluation is ~`Q` Phi evaluations
  per occurrence/background point relative to the Gaussian kernel.
* `likelihood = "skew_t"`: presence-only fit (no `env_m` required).
* `likelihood = "skew_t_weighted"`: paper Eq. 5 with the NCST density
  in place of the Gaussian, plus the same ridge prior on `log sigma`
  used by `"weighted"` and `"skew_normal_weighted"`. Defaults inherit
  from `"skew_normal_weighted"`; pass `prior_log_sigma_lambda = 0` for
  pure paper-faithful skew-t MLE.
* New theta layout: `[mu (p), log_sigma (p), v (p(p-1)/2), alpha (p),
  log_r (1)]` -- one extra parameter relative to the skew-normal
  layout. The Sobol bound for `log r` defaults to
  `[log 2, log 100]` (so `r` in `[2, 100]`), with the Sobol centre at
  `log 10`.
* Both families currently use finite-difference gradients (the C++
  kernel evaluates the GL quadrature on every call; analytic gradients
  through the integral are deferred).
* `geom_nicher_isosuitability()` is fully aware of the skew-t theta
  layout and contours the NCST density on a 2-D grid using a 20-node
  Gauss--Laguerre approximation.

### Skew-normal niche likelihoods (`skew_normal`, `skew_normal_weighted`)

* Two new likelihood families on `optimize_niche()` that fit a
  multivariate skew-normal niche
  (Azzalini & Capitanio 1999, *J. R. Stat. Soc. Ser. B* **61**(3),
  579--602) instead of the symmetric Gaussian:
  ```
  S(x) ∝ φ_p(x − μ; Σ) · Φ( Σ_k α_k · (x_k − μ_k) / σ_k )
  ```
  Adds a length-`p` skewness vector `α ∈ ℝ^p` to the existing
  `(μ, Σ)` parameters; `α_k = 0` for all `k` recovers the Gaussian
  niche exactly. The skew niche is useful when the environmental
  conditions at occurrences are not symmetric around the mode --
  e.g. species that tolerate higher precipitation but avoid drier
  conditions, or vice-versa.
* `likelihood = "skew_normal"`: presence-only fit (no `env_m`
  required). Closed-form Azzalini density; finite-difference
  gradient (analytic Azzalini gradients deferred).
* `likelihood = "skew_normal_weighted"`: paper Eq. 5 with the SN
  density in place of the Gaussian, plus the same ridge prior on
  `log sigma` used by `"weighted"` (default `lambda = 1.0`,
  `prior_log_sigma_center = log(sd(env_occ))`). Defaults inherit
  from `"weighted"`; pass `prior_log_sigma_lambda = 0` for pure
  paper-faithful skew-normal MLE (not recommended on multimodal M).
* The skew families currently require `backend = "cpp"`; the legacy
  R backend has no pure-R reference implementation for them.
* Recommended `α` interpretation: `|α_k| <= 3` covers the practically
  useful skewness range per Azzalini & Capitanio (1999) §5; the
  Sobol-start machinery samples `α_k` in `[-3, 3]`.

### `geom_nicher_isosuitability()`: contour layer for arbitrary niche
geometries

* New `ggplot2` layer `geom_nicher_isosuitability(model, level, ...)`
  that draws iso-suitability contours `S(x) = c` by evaluating the
  fitted suitability function on a 2-D environmental grid and
  contouring the result. Unlike `geom_nicher_ellipse()`, which traces
  analytical ellipses tied to the Gaussian niche geometry, this layer
  works for any likelihood family supported by `optimize_niche()`,
  including the skew-normal families where iso-suitability sets are
  **not** ellipses (the Φ(z) factor breaks ellipse symmetry).
* Suitability is normalised so `S(μ) = 1` at the SN/Gaussian location
  parameter, so `level = 0.5` always means "regions where the niche
  is at least 50 % as suitable as the reference centre".
* For Gaussian fits the contour at level `c` is a closed ellipse
  identical to `geom_nicher_ellipse(level = c, level_type =
  "suitability")`. For skew-normal fits the contours bunch on the
  side that `α` pulls suitability toward and stretch on the opposite
  side -- using `geom_nicher_ellipse()` on a skew-normal fit would
  draw an ellipse that misrepresents the model.

### ggplot2 plotting engine for 2-D fitted niches (E-space)

* New `autoplot.nicher()` S3 method on `ggplot2::autoplot` returns a
  `ggplot` object visualizing a fitted 2-D `nicher` model in
  environmental space. Composes a faint background point cloud
  (`env_m`), an occurrence cloud (`env_occ`), and the model-implied
  iso-suitability ellipses; axes are labelled by `object$var_names`.
* New composable layer constructors — each returns a single self-
  contained `ggplot2::layer` and is stackable with `+` on a fresh
  `ggplot()`:
    * `geom_nicher_background(env_m, ...)` — background environment.
    * `geom_nicher_occ(env_occ, ...)` — occurrence cloud.
    * `geom_nicher_ellipse(model, level, level_type, n, ...)` — niche
      geometry; default `level = c(0.95, 0.5, 0.05)` interpreted as
      iso-suitability contours `S(x) = level` (paper-aligned, Jiménez
      et al. 2022 Eq. 2). Pass `level_type = "chisq"` for the textbook
      χ² confidence-ellipse interpretation.
* New `nicher_compare_plot(models, env_occ, env_m, ...)` convenience
  wrapper for the two multi-model contexts called out in the design
  spec:
    * Single species, multiple models — `env_occ` is shared.
    * Multiple species, multiple models — `env_occ` is a named list
      whose names match `names(models)`; ellipse colour identifies the
      model and occurrences are coloured to match.
* Plotting is restricted to 2-D models (`length(object$var_names) ==
  2`). Higher-dimensional fits are rejected with a clear error; no
  projection or marginalization is performed.
* Fitted `nicher` objects do not store environmental data — `env_occ`
  and `env_m` must be supplied explicitly to every plotting call.
* `ggplot2 (>= 3.0.0)` is added as a *Suggests* dependency. nicher
  installs and runs without ggplot2; the plotting methods register at
  load time only when ggplot2 is available, and each entry point
  guards with a clear error message otherwise. The `linewidth` /
  `size` argument rename across ggplot2 3.4.0 is handled
  transparently.

# nicher 2.2.3

## Bug fixes (Windows, follow-up)

* The 2.2.2 attempt at forcing `RcppParallel`'s header-only TinyThread
  backend was incomplete: `src/Makevars{,.win}` set `PKG_CXXFLAGS` to
  ```
  -DRCPP_PARALLEL_USE_TBB=0 $(shell ... RcppParallel::CxxFlags())
  ```
  but on Windows `RcppParallel::CxxFlags()` expands to
  `-I<...>/RcppParallel/include -DRCPP_PARALLEL_USE_TBB=1`, so g++ saw
  ```
  -DRCPP_PARALLEL_USE_TBB=0 ... -DRCPP_PARALLEL_USE_TBB=1
  ```
  with the *last* `-D` winning, silently turning TBB back on at the
  source level. The compiled DLL then failed to load with the original
  `LoadLibrary failure: The specified module could not be found`
  because tbb.dll / tbbmalloc.dll live inside the `RcppParallel`
  install dir and are not on the default DLL search path.
* Both `Makevars` files now set `PKG_CXXFLAGS = -DRCPP_PARALLEL_USE_TBB=0`
  exactly once. `RcppParallel`'s include directory is supplied
  automatically by `LinkingTo: RcppParallel` (declared in DESCRIPTION),
  so dropping the `RcppParallel::CxxFlags()` `$(shell ...)` substitution
  loses nothing.
* Verified locally on Linux: the per-source compile line now contains a
  single `-DRCPP_PARALLEL_USE_TBB=0`, no `<command-line>: warning:
  "RCPP_PARALLEL_USE_TBB" redefined`, and `library(nicher)` loads
  cleanly. Windows users should see the same.

# nicher 2.2.2

## Bug fixes (Windows)

* `nicher.dll` no longer fails to load on Windows with
  `LoadLibrary failure: The specified module could not be found`.
  Both `src/Makevars` and `src/Makevars.win` now drop
  `RcppParallel::RcppParallelLibs()` from `PKG_LIBS`. The TinyThread
  backend (`-DRCPP_PARALLEL_USE_TBB=0`) was already enabled at compile
  time but the linker line still referenced TBB / `tbbmalloc`, whose
  DLLs live inside the `RcppParallel` install directory and are not on
  Windows' default DLL search path. Forcing the header-only TinyThread
  backend at both compile *and* link time eliminates the runtime
  dependency on `tbb.dll` / `tbbmalloc.dll`. Same change is applied to
  `Makevars` (Linux/macOS) for consistency.

## Validation & UX

* `optimize_niche()` now validates `eta` upfront, rejecting non-numeric,
  non-positive, non-finite, or non-scalar values *before* any compute is
  spent on a multi-start optimization. Previously the check happened in
  `new_nicher()` after the optimization had already run.
* `print.nicher()` now surfaces the `eta` value used at fit time and the
  variable names stored on the object (when present).

## Dependency hygiene

* `RcppEigen` removed from the `Imports:` field — it is LinkingTo only,
  with no R-level `::` usage. Clears the
  "Namespaces in Imports not imported from: RcppEigen" NOTE.
* `RcppParallel` retained in `Imports:` (it is called at runtime as
  `RcppParallel::defaultNumThreads()`) but now declared in `NAMESPACE`
  via `@importFrom RcppParallel defaultNumThreads` on
  `habitat_suitability()`. Clears the matching NOTE for `RcppParallel`.

## Documentation

* `R CMD check`: 5 WARNINGs cleared, down from 5 W + 3 N to 0 W + 2 N
  (installed size and `SystemRequirements: GNU make`, both unavoidable
  and informational).
* `R/benchmark_optimize_niche.R`: replaced unknown `\lifecycle{...}`
  Rd macro with the canonical
  `\Sexpr[results=rd, stage=render]{lifecycle::badge("...")}` recipe.
* `R/optimize_niche.R`: replaced bare `\cdot` with `*` inside `\code{}`
  (Rd does not honour `\cdot` outside `\eqn{}`).
* `R/niche_kde_bias_corrected.R`: replaced dangling `\link{create_niche_obj_ptr}`
  with plain `\code{}` — the helper is internal and has no Rd file.
* `vignettes/nicher-intro.Rmd` now builds correctly during
  `R CMD build`, populating `inst/doc/`.

# nicher 2.2.1

## Bug fixes

* `predict.nicher()` now reconstructs the correlation matrix using the same
  `eta` value that `optimize_niche()` used at fit time, rather than always
  passing `eta = 1` to `cvine_cholesky()`. `eta` is persisted on the
  returned `nicher` object as a new field. Legacy `nicher` objects produced
  before 2.2.1 (with no `eta` field) continue to predict against the
  default `eta = 1`. Reported by Devin Review on PR #34.

## Test-suite fixes

* `test-benchmark-optimizers.R`: the 3D R-vs-C++ test now verifies kernel
  parity (`fn_r(theta) == fn_cpp(theta)` at each backend's optimum) rather
  than asserting the two optimizers find the same local minimum, which is
  not guaranteed for a non-convex objective when the FD gradient differs
  between the interpreter and the compiled kernel.
* `test-benchmark-optimizers.R`: the XPtr-backend convergence test now
  accepts ucminf code 4 ("zero step from line search") as a valid
  termination, matching the convention used by the 2D R-vs-C++ test.
* `test-benchmark-optimizers.R`: removed `label = ` argument from
  `expect_s3_class()` (unsupported in installed testthat).
* `niche_kde_bias_corrected()` and `niche_presence_only()` now emit informative
  errors for `precomp_w_den` length mismatches and non-finite `start`
  values; the corresponding `test-niche-wrappers.R` cases now exercise
  these paths via named arguments rather than relying on positional
  argument order.

# nicher 2.2.0

## New features

* New `habitat_suitability(param, env, ...)` function: evaluates the
  standardized multivariate-normal suitability map of Jimenez et al.
  (2022, Eq. 2) over a multi-layer
  [`terra::SpatRaster`][terra::SpatRaster] using a streaming
  [`RcppParallel`][RcppParallel] kernel. Memory is bounded by the size
  of one raster block (continental rasters never materialise in R);
  `NA` cells are masked, compacted before the C++ kernel, and
  scattered back into the output. Forwards `output`, `overwrite`, and
  `wopt` to [`terra::writeStart()`][terra::writeStart].
* New `predict.nicher()` S3 method: turns a fitted `nicher` object
  from `optimize_niche()` into a habitat-suitability raster.
  Reconstructs `(mu, Sigma)` from `best$theta` via `cvine_cholesky()`
  and dispatches to `habitat_suitability()`. Reorders the input
  raster by layer name to match the variable order used at fit
  time — supports calling `predict()` on future-climate stacks
  whose layers may be in any order.
* `new_nicher()` and `optimize_niche()` now record the column names
  of `env_occ` on the returned object as `var_names`, used by
  `predict.nicher()` for name-based layer matching.
* New low-level entry point `niche_suitability_cpp()`: parallel
  triangular-solve kernel that takes a flat column-major
  environmental buffer and a precomputed `L_inv` (lower-Cholesky
  inverse of `Sigma`) and returns suitability values pixel-by-pixel.

## Dependencies

* Adds `terra`, `RcppParallel`, and `checkmate` to `Imports`; adds
  `RcppParallel` to `LinkingTo`.

# nicher 2.1.0

## New features

* `optimize_niche()` gains a `backend` argument:
  - `"cpp"` (default): optimizes via `ucminfcpp::ucminf_xptr()` with a pure-C++
    objective. The full theta-unpacking, C-vine Cholesky, log-likelihood, and
    gradient are evaluated in compiled C++ with no R-callback overhead.
  - `"r"`: optimizes via `ucminf::ucminf()` with the legacy R-level objective
    functions. Provided ONLY for side-by-side benchmarking and emits a
    `lifecycle::deprecate_soft()` warning. Will be removed in nicher 2.2.0.
* `optimize_niche()` gains a `grad` argument selecting the gradient strategy:
  `"auto"` (default), `"analytic"`, `"central"`, or `"forward"`. With
  `backend = "cpp"` and the weighted likelihood, `"analytic"` uses a hybrid
  analytic / FD gradient: closed-form derivatives over the location (mu) and
  log-scale (log sigma) blocks, central finite differences over the C-vine
  partial-correlation block.
* `optimize_niche()` gains `m_subsample`, `m_kde_subsample`, and `seed`
  arguments controlling the KDE sampling policy for the weighted model.
  Defaults cap both subsamples at 10 000 distinct background combinations
  (the maximum the model can usefully exploit) and warn below the heuristic
  representative-sample floor `max(500, 50 * 2^p)`.
* New `benchmark_optimize_niche()` function runs both backends from identical
  Sobol starts and returns timing / convergence statistics. Will be removed
  alongside `backend = "r"` in 2.2.0.

## Internal

* New pure-C++ math-scale kernels:
  - `loglik_niche_math_presence_only_cpp()`
  - `loglik_niche_math_kde_bias_corrected_cpp()` (with precomputed KDE weights)
  - `loglik_niche_math_kde_bias_corrected_grad_cpp()` (hybrid analytic gradient)
  These power the C++ backend and bypass `Rcpp::NumericMatrix` allocations
  inside the optimizer's inner loop.
* New Eigen-native `nicher::cvine_cholesky_eigen()` helper.
* `create_niche_obj_ptr()` accepts `precomp_w_occ` and `grad` arguments and
  wraps the closure body in `try / catch` so that a failed evaluation
  surfaces as `f = +Inf, g = 0` rather than aborting R.

# nicher 2.0.0

* Initial public release.
