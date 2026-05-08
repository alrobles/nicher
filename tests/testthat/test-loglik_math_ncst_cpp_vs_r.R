# Stage 3 QA gate: cross-language parity between the pure-R reference
# implementation (loglik_niche_math_ncst_r) and the C++/Eigen kernel
# (loglik_niche_math_ncst_cpp) for the Non-Central Skew t (NCST) density.

test_that("NCST presence-only: C++ matches R reference (p = 2)", {

  set.seed(42)
  p <- 2
  n <- 30
  env_occ <- matrix(rnorm(n * p), ncol = p)
  colnames(env_occ) <- paste0("var", seq_len(p))

  n_v <- p * (p - 1) / 2
  theta <- c(
    runif(p, -1, 1),         # xi
    runif(p, -0.5, 0.5),     # log_sigma
    runif(n_v, -1, 1),       # v (C-vine partials)
    runif(p, -2, 2),         # alpha
    log(10)
  )

  res_r <- loglik_niche_math_ncst_r(theta, env_occ, eta = 1.0)
  res_cpp <- loglik_niche_math_ncst_cpp(theta, env_occ, eta = 1.0)

  expect_equal(res_cpp, res_r, tolerance = 1e-10)
})

test_that("NCST presence-only: C++ matches R reference (p = 3)", {
  set.seed(123)
  p <- 3
  n <- 25
  env_occ <- matrix(rnorm(n * p, sd = 2), ncol = p)
  colnames(env_occ) <- paste0("var", seq_len(p))

  n_v <- p * (p - 1) / 2
  theta <- c(
    c(0.5, -0.3, 1.0),      # xi
    c(0.2, -0.1, 0.3),      # log_sigma
    runif(n_v, -0.5, 0.5),  # v
    c(1.5, -1.0, 0.5),      # alpha
    log(5)
  )

  res_r <- loglik_niche_math_ncst_r(theta, env_occ, eta = 1.0)
  res_cpp <- loglik_niche_math_ncst_cpp(theta, env_occ, eta = 1.0)

  expect_equal(res_cpp, res_r, tolerance = 1e-10)
})

test_that("NCST presence-only: multiple parameter sets yield parity", {
  set.seed(999)
  p <- 2
  n <- 20
  env_occ <- matrix(rnorm(n * p), ncol = p)
  n_v <- p * (p - 1) / 2

  for (trial in seq_len(5)) {
    theta <- c(
      rnorm(p),                # xi
      rnorm(p, sd = 0.5),     # log_sigma
      runif(n_v, -2, 2),      # v
      rnorm(p, sd = 1.5),     # alpha
      log(runif(1, 2, 50))    # log_r
    )
    res_r <- loglik_niche_math_ncst_r(theta, env_occ, eta = 1.0)
    res_cpp <- loglik_niche_math_ncst_cpp(theta, env_occ, eta = 1.0)
    expect_equal(res_cpp, res_r, tolerance = 1e-10,
                 label = paste("trial", trial))
  }
})

test_that("NCST presence-only: extreme r values yield parity", {
  set.seed(7)
  p <- 2
  n <- 15
  env_occ <- matrix(rnorm(n * p), ncol = p)
  n_v <- p * (p - 1) / 2

  for (log_r_val in c(log(2), log(3), log(50), log(100))) {
    theta <- c(
      c(0, 0),               # xi
      c(0, 0),
      0,
      c(1, -1),              # alpha
      log_r_val              # log_r
    )
    res_r <- loglik_niche_math_ncst_r(theta, env_occ, eta = 1.0)
    res_cpp <- loglik_niche_math_ncst_cpp(theta, env_occ, eta = 1.0)
    expect_equal(res_cpp, res_r, tolerance = 1e-10,
                 label = paste("log_r =", round(log_r_val, 2)))
  }
})

test_that("NCST presence-only: zero alpha yields symmetric (no skew)", {
  set.seed(11)
  p <- 2
  n <- 20
  env_occ <- matrix(rnorm(n * p), ncol = p)
  n_v <- p * (p - 1) / 2

  theta <- c(
    c(0, 0),         # xi
    c(0, 0),         # log_sigma
    0,               # v
    c(0, 0),         # alpha = 0 (no skew)
    log(10)          # log_r
  )

  res_r <- loglik_niche_math_ncst_r(theta, env_occ, eta = 1.0)
  res_cpp <- loglik_niche_math_ncst_cpp(theta, env_occ, eta = 1.0)
  expect_equal(res_cpp, res_r, tolerance = 1e-10)
  expect_true(is.finite(res_cpp))
})

test_that("NCST weighted: C++ matches R for basic case", {
  skip("Weighted R reference not yet implemented")
})

test_that("NCST returns finite for typical ecological data", {
  set.seed(55)
  p <- 2
  n <- 50
  env_occ <- matrix(c(rnorm(n, mean = 20, sd = 5),
                       rnorm(n, mean = 500, sd = 100)),
                     ncol = p)
  colnames(env_occ) <- c("bio1", "bio12")
  n_v <- p * (p - 1) / 2

  theta <- c(
    c(20, 500),              # xi near data means
    c(log(5), log(100)),     # log_sigma near data SDs
    0,                       # v
    c(0.5, -0.5),           # moderate skew
    log(10)
  )

  res_r <- loglik_niche_math_ncst_r(theta, env_occ, eta = 1.0)
  res_cpp <- loglik_niche_math_ncst_cpp(theta, env_occ, eta = 1.0)
  expect_equal(res_cpp, res_r, tolerance = 1e-10)
  expect_true(is.finite(res_r))
  expect_true(is.finite(res_cpp))
})
