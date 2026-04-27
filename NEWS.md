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
