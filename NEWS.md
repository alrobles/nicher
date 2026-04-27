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
* `R/niche_weighted.R`: replaced dangling `\link{create_niche_obj_ptr}`
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
* `niche_weighted()` and `niche_presence_only()` now emit informative
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
  - `loglik_niche_math_weighted_cpp()` (with precomputed KDE weights)
  - `loglik_niche_math_weighted_grad_cpp()` (hybrid analytic gradient)
  These power the C++ backend and bypass `Rcpp::NumericMatrix` allocations
  inside the optimizer's inner loop.
* New Eigen-native `nicher::cvine_cholesky_eigen()` helper.
* `create_niche_obj_ptr()` accepts `precomp_w_occ` and `grad` arguments and
  wraps the closure body in `try / catch` so that a failed evaluation
  surfaces as `f = +Inf, g = 0` rather than aborting R.

# nicher 2.0.0

* Initial public release.
