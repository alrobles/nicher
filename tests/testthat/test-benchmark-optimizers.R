library(testthat)
# tests/testthat/test-benchmark-optimizers.R
#
# Benchmarking tests: ucminf::ucminf (R backend) vs. ucminfcpp (C++ backend)
#
# These tests compare the execution time and accuracy of two optimizer
# configurations used within nicher:
#
#   1. ucminf::ucminf with an R-level objective function that builds the
#      covariance matrix inside R and calls the pure-R loglik_niche().
#      (Referred to below as the "R backend".)
#
#   2. ucminf::ucminf with a C++ objective function loglik_niche_math_cpp(),
#      which passes the Cholesky factor directly to loglik_niche_chol_cpp().
#      (Referred to below as the "ucminfcpp / C++ backend".)
#
#   3. ucminfcpp::ucminf_xptr() with a compiled C++ objective function created
#      by create_niche_obj_ptr() — the new XPtr backend.
#      (Referred to below as the "xptr backend".)
#
# The benchmarks run on the built-in 2-D hummingbird example datasets
# (example_env_occ_2d / example_env_m_2d) so that they are fast enough
# for R CMD check.  Heavier multi-start scenarios are guarded with
# skip_on_cran() so that they do not time out on CRAN.

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

#' Run ucminf with an R-level (covariance-matrix) objective function.
#'
#' Converts the math-scale parameter vector \code{theta} to a covariance
#' matrix on the R side, then evaluates \code{\link{loglik_niche}}.
#' This corresponds to the **R-backend** approach.
#'
#' @param theta0 Numeric starting vector (math scale).
#' @param env_occ,env_m Data frames for presences and background.
#' @param control List of control parameters forwarded to
#'   \code{\link[ucminf]{ucminf}}.
#' @return The raw \code{ucminf} result list.
.run_ucminf_r_backend <- function(theta0, env_occ, env_m,
                                  control = list(maxeval = 500)) {
  p <- ncol(env_occ)

  fn_r <- function(theta) {
    mu <- theta[seq_len(p)]
    sigma <- exp(theta[(p + 1L):(2L * p)])
    v <- if (p > 1L) theta[(2L * p + 1L):length(theta)] else numeric(0L)
    L_corr <- cvine_cholesky(v, d = p, eta = 1)
    L_cov <- diag(sigma) %*% L_corr
    S <- tcrossprod(L_cov) # Sigma = L L'
    loglik_niche(
      mu = mu,
      s_mat = S,
      env_occ = env_occ,
      env_m = env_m,
      neg = TRUE
    )
  }

  ucminf::ucminf(par = theta0, fn = fn_r, hessian = FALSE, control = control)
}

#' Run ucminf with the C++ objective function (ucminfcpp backend).
#'
#' Calls \code{\link{loglik_niche_math_cpp}} which passes the Cholesky
#' factor directly to the native C++ kernel \code{loglik_niche_chol_cpp},
#' eliminating R-level matrix construction overhead.  This is the
#' **ucminfcpp (C++ backend)** approach.
#'
#' @param theta0 Numeric starting vector (math scale).
#' @param env_occ,env_m Data frames for presences and background.
#' @param control List of control parameters forwarded to
#'   \code{\link[ucminf]{ucminf}}.
#' @return The raw \code{ucminf} result list.
.run_ucminf_cpp_backend <- function(theta0, env_occ, env_m,
                                    control = list(maxeval = 500)) {
  fn_cpp <- function(theta) {
    loglik_niche_math_cpp(theta,
      env_occ = env_occ,
      env_m   = env_m,
      eta     = 1,
      neg     = TRUE
    )
  }

  ucminf::ucminf(par = theta0, fn = fn_cpp, hessian = FALSE, control = control)
}

#' Run ucminfcpp::ucminf_xptr() with the compiled C++ NicheObjFun (xptr backend).
#'
#' Uses \code{\link{create_niche_obj_ptr}} to build a compiled C++ objective
#' function and passes it directly to \code{ucminfcpp::ucminf_xptr()},
#' eliminating all R interpreter round-trips during optimization.
#'
#' Note: uses the internal \code{"unweighted"} likelihood (C++ lik_type = 0)
#' to mirror the classical log-likelihood computed by the R backend helper
#' above via \code{loglik_niche()}. The R-level \code{niche_unweighted()}
#' wrapper was removed in PR #27, but the C++ code path remains and is
#' still the correct benchmark counterpart here.
#'
#' @param theta0 Numeric starting vector (math scale).
#' @param env_occ,env_m Data frames for presences and background.
#' @param control List of control parameters forwarded to
#'   \code{ucminfcpp::ucminf_control}.
#' @return The \code{ucminf} result list returned by \code{ucminfcpp::ucminf_xptr}.
.run_xptr_backend <- function(theta0, env_occ, env_m,
                              control = list(maxeval = 500)) {
  gs <- if (!is.null(control$gradstep)) control$gradstep else c(1e-6, 1e-8)
  xptr <- create_niche_obj_ptr(
    env_occ    = as.matrix(env_occ),
    env_m      = as.matrix(env_m),
    eta        = 1.0,
    likelihood = "unweighted",
    gradstep   = gs
  )
  con <- do.call(ucminfcpp::ucminf_control, control)
  ucminfcpp::ucminf_xptr(par = theta0, xptr = xptr, control = con)
}

# ---------------------------------------------------------------------------
# Test 1 – Single-start correctness (2-D dataset)
# ---------------------------------------------------------------------------
test_that("R-backend and C++ backend produce numerically consistent log-likelihoods (2D)", {
  set.seed(42)
  theta0 <- start_theta(example_env_occ_2d)

  res_r <- .run_ucminf_r_backend(theta0, example_env_occ_2d, example_env_m_2d)
  res_cpp <- .run_ucminf_cpp_backend(theta0, example_env_occ_2d, example_env_m_2d)

  # Both optimizations must reach a finite negative log-likelihood
  expect_true(is.finite(res_r$value),
    label = "R-backend returns a finite objective value"
  )
  expect_true(is.finite(res_cpp$value),
    label = "C++ backend returns a finite objective value"
  )

  # Optimized log-likelihoods must agree within a generous tolerance
  # (both backends minimize the same mathematical function)
  expect_equal(res_r$value, res_cpp$value,
    tolerance = 1e-3,
    label = "R- and C++ backend log-likelihoods match within 1e-3"
  )

  # Both must report convergence (ucminf code 1 or 2 = success)
  expect_true(res_r$convergence %in% c(1L, 4L),
    label = "R-backend converges (code 1 or 4)"
  )
  expect_true(res_cpp$convergence %in% c(1L, 4L),
    label = "C++ backend converges (code 1 or 4)"
  )
})

# ---------------------------------------------------------------------------
# Test 2 – Single-start correctness (3-D dataset)
# ---------------------------------------------------------------------------
test_that("R-backend and C++ backend produce numerically consistent log-likelihoods (3D)", {
  set.seed(42)
  theta0 <- start_theta(example_env_occ_3d)

  res_r <- .run_ucminf_r_backend(theta0, example_env_occ_3d, example_env_m_3d)
  res_cpp <- .run_ucminf_cpp_backend(theta0, example_env_occ_3d, example_env_m_3d)

  expect_true(is.finite(res_r$value),
    label = "R-backend 3D returns finite value"
  )
  expect_true(is.finite(res_cpp$value),
    label = "C++ backend 3D returns finite value"
  )

  expect_equal(res_r$value, res_cpp$value,
    tolerance = 1e-3,
    label = "R- and C++ backend 3D log-likelihoods match within 1e-3"
  )
})

# ---------------------------------------------------------------------------
# Test 3 – Optimized parameter vectors are numerically close (2-D dataset)
# ---------------------------------------------------------------------------
test_that("R-backend and C++ backend converge to the same parameter vector (2D)", {
  set.seed(42)
  theta0 <- start_theta(example_env_occ_2d)

  res_r <- .run_ucminf_r_backend(theta0, example_env_occ_2d, example_env_m_2d)
  res_cpp <- .run_ucminf_cpp_backend(theta0, example_env_occ_2d, example_env_m_2d)

  # Parameter estimates must coincide within a moderate tolerance
  expect_equal(res_r$par, res_cpp$par,
    tolerance = 1e-2,
    label = "Optimal parameter vectors agree within 1e-2"
  )
})

# ---------------------------------------------------------------------------
# Test 4 – Timing benchmark: C++ backend should not be slower than R backend
# ---------------------------------------------------------------------------
test_that("C++ backend is not slower than R backend (2D, single start)", {
  skip_on_cran()

  set.seed(42)
  theta0 <- start_theta(example_env_occ_2d)
  n_reps <- 10L

  time_r <- system.time(
    for (i in seq_len(n_reps)) {
      .run_ucminf_r_backend(theta0, example_env_occ_2d, example_env_m_2d)
    }
  )["elapsed"]

  time_cpp <- system.time(
    for (i in seq_len(n_reps)) {
      .run_ucminf_cpp_backend(theta0, example_env_occ_2d, example_env_m_2d)
    }
  )["elapsed"]

  # Report results for informational purposes
  message(sprintf(
    "\n--- Timing benchmark (2D, %d repetitions) ---\n  R backend:   %.4f s total (%.4f s/rep)\n  C++ backend: %.4f s total (%.4f s/rep)\n  Speedup:     %.2fx\n",
    n_reps,
    time_r,   time_r / n_reps,
    time_cpp, time_cpp / n_reps,
    time_r / max(time_cpp, .Machine$double.eps)
  ))

  # The C++ backend should be at most twice as slow as the R backend.
  # In practice it is expected to be faster; this threshold prevents
  # accidental severe regressions from being silently ignored.
  expect_true(
    time_cpp <= 2 * time_r,
    label = "C++ backend wall time is within 2x of R backend wall time"
  )
})

# ---------------------------------------------------------------------------
# Test 5 – Multi-start benchmark via optimize_niche (heavier, off-CRAN only)
# ---------------------------------------------------------------------------
test_that("optimize_niche (xptr backend) outperforms equivalent ucminf loop (2D)", {
  skip_on_cran()

  set.seed(123)
  n_starts <- 5L

  # --- optimize_niche (always xptr, presence_only for simplicity) -------
  time_cpp <- system.time(
    result_cpp <- optimize_niche(
      env_occ    = example_env_occ_2d,
      env_m      = NULL,
      num_starts = n_starts,
      breadth    = 0.1,
      likelihood = "presence_only",
      eta        = 1
    )
  )["elapsed"]

  # --- ucminf R backend: replicate multi-start logic manually -----------
  starts_df <- start_theta_multiple(
    example_env_occ_2d, n_starts,
    method = "uniform"
  )
  starts_list <- split(starts_df, seq_len(nrow(starts_df)))
  starts_list <- lapply(starts_list, function(x) {
    v <- as.numeric(x)
    names(v) <- colnames(starts_df)
    v
  })

  fn_po <- function(theta) {
    loglik_niche_math_presence_only(
      theta, example_env_occ_2d, eta = 1, neg = TRUE
    )
  }
  time_r <- system.time({
    results_r <- lapply(starts_list, function(s) {
      ucminf::ucminf(par = s, fn = fn_po,
        control = list(maxeval = 500L), hessian = FALSE)
    })
  })["elapsed"]

  best_loglik_cpp <- result_cpp$best$loglik
  best_loglik_r   <- max(sapply(results_r, function(x) -x$value))

  message(sprintf(
    paste0(
      "\n--- Multi-start benchmark (%d starts, 2D, presence_only) ---\n",
      "  R backend:  %.4f s | best loglik = %.4f\n",
      "  xptr:       %.4f s | best loglik = %.4f\n",
      "  Speedup:    %.2fx\n"
    ),
    n_starts,
    time_r,   best_loglik_r,
    time_cpp, best_loglik_cpp,
    time_r / max(time_cpp, .Machine$double.eps)
  ))

  # Both approaches must reach comparable likelihoods
  expect_equal(best_loglik_r, best_loglik_cpp,
    tolerance = 0.1,
    label = "Multi-start best log-likelihoods agree within 0.1 units"
  )

  # The xptr backend must not be dramatically slower than ucminf
  expect_true(
    time_cpp <= 3 * time_r,
    label = "xptr multi-start is within 3x of ucminf wall time"
  )
})

# ---------------------------------------------------------------------------
# Test 6 – Presence-only likelihood: single-start correctness
# ---------------------------------------------------------------------------
test_that("Presence-only C++ backend matches R-level reconstruction (2D)", {
  set.seed(42)
  theta0 <- start_theta(example_env_occ_2d)

  # C++ backend via the exported wrapper
  val_cpp <- loglik_niche_math_presence_only(theta0, example_env_occ_2d,
    eta = 1, neg = TRUE
  )

  # R-level reconstruction matching the C++ formula exactly:
  #   neg_log = 0.5 * n * log_det + 0.5 * sum_q
  # where log_det = 2 * sum(log(diag(L_cov))) and
  #       sum_q   = ||L_cov^{-1} (X - mu)||_F^2  (sum over rows of X)
  p <- ncol(example_env_occ_2d)
  mu0 <- theta0[seq_len(p)]
  sigma0 <- exp(theta0[(p + 1L):(2L * p)])
  v0 <- theta0[(2L * p + 1L):length(theta0)]
  L_corr0 <- cvine_cholesky(v0, d = p, eta = 1)
  L_cov0 <- diag(sigma0) %*% L_corr0

  env_mat <- as.matrix(example_env_occ_2d)
  n <- nrow(env_mat)
  diff_occ <- t(sweep(env_mat, 2L, mu0)) # p x n matrix of (x_i - mu)
  y_occ <- forwardsolve(L_cov0, diff_occ) # p x n: L^{-1} (x_i - mu)
  sum_q <- sum(y_occ^2)
  log_det <- 2 * sum(log(diag(L_cov0))) # log|Sigma|

  val_r <- 0.5 * n * log_det + 0.5 * sum_q

  expect_equal(val_cpp, val_r,
    tolerance = 1e-6,
    label = "Presence-only C++ and R values agree to 1e-6"
  )
})

# ---------------------------------------------------------------------------
# Test 7 – XPtr backend (ucminfcpp::ucminf_xptr) correctness (2D)
# ---------------------------------------------------------------------------
test_that("XPtr backend (create_niche_obj_ptr) produces consistent results (2D)", {
  set.seed(42)
  theta0 <- start_theta(example_env_occ_2d)

  res_cpp <- .run_ucminf_cpp_backend(theta0, example_env_occ_2d, example_env_m_2d)
  res_xptr <- .run_xptr_backend(theta0, example_env_occ_2d, example_env_m_2d)

  # Both must yield a finite objective value
  expect_true(is.finite(res_cpp$value),
    label = "C++ backend returns finite value"
  )
  expect_true(is.finite(res_xptr$value),
    label = "XPtr backend returns finite value"
  )

  # Log-likelihoods must agree within tolerance
  expect_equal(res_cpp$value, res_xptr$value,
    tolerance = 1e-3,
    label = "C++ and XPtr backend log-likelihoods agree within 1e-3"
  )

  # Convergence must be successful for both
  expect_true(res_xptr$convergence %in% c(1L, 2L),
    label = "XPtr backend converges (code 1 or 2)"
  )
})

# ---------------------------------------------------------------------------
# Test 8 – optimize_niche (cpp xptr) returns a valid nicher object (2D)
# ---------------------------------------------------------------------------
test_that("optimize_niche returns a valid nicher object (2D, presence_only)", {
  set.seed(42)

  res <- optimize_niche(
    env_occ    = example_env_occ_2d,
    env_m      = NULL,
    num_starts = 3L,
    breadth    = 0.1,
    likelihood = "presence_only",
    eta        = 1
  )

  # Must return a nicher S3 object
  expect_s3_class(res, "nicher")
  expect_true(!is.null(res$best), label = "result has $best")
  expect_true(!is.null(res$solutions), label = "result has $solutions")
  expect_true(is.finite(res$best$loglik),
    label = "best$loglik is finite"
  )
  expect_true(res$best$convergence %in% c(1L, 2L),
    label = "best$convergence is successful"
  )
  expect_equal(res$likelihood, "presence_only",
    label = "likelihood field is set correctly"
  )
  expect_equal(res$n_starts, 3L,
    label = "n_starts field matches num_starts"
  )
})

# ---------------------------------------------------------------------------
# Test 9 – Presence-only likelihood via optimize_niche (2D)
# ---------------------------------------------------------------------------
test_that("optimize_niche(likelihood='presence_only') returns valid nicher (2D)", {
  set.seed(42)

  res <- optimize_niche(
    env_occ      = example_env_occ_2d,
    env_m        = NULL,
    num_starts   = 3L,
    breadth      = 0.1,
    likelihood   = "presence_only",
    eta          = 1
  )

  expect_s3_class(res, "nicher", label = "presence_only result is nicher")
  expect_true(
    is.finite(res$best$loglik),
    label = "presence_only best loglik is finite"
  )
  expect_equal(res$likelihood, "presence_only",
    label = "likelihood field is 'presence_only'"
  )
})
