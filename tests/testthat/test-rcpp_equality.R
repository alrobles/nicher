# Tests: R vs C++ log-likelihood equality
# Validates that loglik_unweighted_math == loglik_unweighted_cpp and
# loglik_presenceonly_math == loglik_presenceonly_cpp on the same inputs.

# ── helpers ─────────────────────────────────────────────────────────────────

make_par <- function(data) {
  par <- get_ellip_par(data)
  par$L <- t(chol(par$S))
  par
}

# ── unweighted model: R vs C++ ───────────────────────────────────────────────

test_that("loglik_unweighted_math equals loglik_unweighted_cpp (package data)", {
  par <- make_par(spOccPnts)
  r_val   <- loglik_unweighted_math(spOccPnts, samMPts, par$mu, par$L)
  cpp_val <- loglik_unweighted_cpp(
    as.matrix(spOccPnts), as.matrix(samMPts), par$mu, par$L
  )
  expect_equal(r_val, cpp_val, tolerance = 1e-8)
})

test_that("loglik_unweighted_math equals loglik_unweighted_cpp (synthetic 2-D data)", {
  set.seed(42)
  occ <- matrix(rnorm(40, mean = c(1, 2), sd = c(0.5, 0.5)), ncol = 2,
                byrow = TRUE)
  bg  <- matrix(rnorm(200, mean = 0, sd = 2), ncol = 2, byrow = TRUE)
  mu  <- colMeans(occ)
  S   <- cov(occ)
  L   <- t(chol(S))

  r_val   <- loglik_unweighted_math(occ, bg, mu, L)
  cpp_val <- loglik_unweighted_cpp(occ, bg, as.vector(mu), L)
  expect_equal(r_val, cpp_val, tolerance = 1e-8)
})

test_that("loglik_unweighted_math equals loglik_unweighted_cpp (3-D data)", {
  set.seed(7)
  occ <- matrix(rnorm(60, mean = c(0, 1, -1), sd = 0.4), ncol = 3,
                byrow = TRUE)
  bg  <- matrix(rnorm(300, mean = 0, sd = 2), ncol = 3, byrow = TRUE)
  mu  <- colMeans(occ)
  S   <- cov(occ)
  L   <- t(chol(S))

  r_val   <- loglik_unweighted_math(occ, bg, mu, L)
  cpp_val <- loglik_unweighted_cpp(occ, bg, as.vector(mu), L)
  expect_equal(r_val, cpp_val, tolerance = 1e-8)
})

test_that("loglik_unweighted R and C++ both return a single finite number", {
  par <- make_par(spOccPnts)
  r_val   <- loglik_unweighted_math(spOccPnts, samMPts, par$mu, par$L)
  cpp_val <- loglik_unweighted_cpp(
    as.matrix(spOccPnts), as.matrix(samMPts), par$mu, par$L
  )
  expect_length(r_val,   1L)
  expect_length(cpp_val, 1L)
  expect_true(is.finite(r_val))
  expect_true(is.finite(cpp_val))
})

# ── presence-only model: R vs C++ ────────────────────────────────────────────

test_that("loglik_presenceonly_math equals loglik_presenceonly_cpp (package data)", {
  par     <- make_par(spOccPnts)
  r_val   <- loglik_presenceonly_math(spOccPnts, samMPts, par$mu, par$S)
  cpp_val <- loglik_presenceonly_cpp(
    as.matrix(spOccPnts), as.matrix(samMPts), par$mu, par$S
  )
  expect_equal(r_val, cpp_val, tolerance = 1e-8)
})

test_that("loglik_presenceonly_math equals loglik_presenceonly_cpp (synthetic 2-D data)", {
  set.seed(99)
  occ <- matrix(rnorm(30, mean = c(2, -1), sd = 0.6), ncol = 2, byrow = TRUE)
  bg  <- matrix(rnorm(150, mean = 0, sd = 3), ncol = 2, byrow = TRUE)
  mu  <- colMeans(occ)
  S   <- cov(occ)

  r_val   <- loglik_presenceonly_math(occ, bg, mu, S)
  cpp_val <- loglik_presenceonly_cpp(occ, bg, as.vector(mu), S)
  expect_equal(r_val, cpp_val, tolerance = 1e-8)
})

test_that("loglik_presenceonly_math equals loglik_presenceonly_cpp (3-D data)", {
  set.seed(13)
  occ <- matrix(rnorm(45, mean = c(1, -2, 0.5), sd = 0.3), ncol = 3,
                byrow = TRUE)
  bg  <- matrix(rnorm(300, mean = 0, sd = 2), ncol = 3, byrow = TRUE)
  mu  <- colMeans(occ)
  S   <- cov(occ)

  r_val   <- loglik_presenceonly_math(occ, bg, mu, S)
  cpp_val <- loglik_presenceonly_cpp(occ, bg, as.vector(mu), S)
  expect_equal(r_val, cpp_val, tolerance = 1e-8)
})

test_that("loglik_presenceonly R and C++ both return a single finite number", {
  par     <- make_par(spOccPnts)
  r_val   <- loglik_presenceonly_math(spOccPnts, samMPts, par$mu, par$S)
  cpp_val <- loglik_presenceonly_cpp(
    as.matrix(spOccPnts), as.matrix(samMPts), par$mu, par$S
  )
  expect_length(r_val,   1L)
  expect_length(cpp_val, 1L)
  expect_true(is.finite(r_val))
  expect_true(is.finite(cpp_val))
})

# ── cross-model sign check ────────────────────────────────────────────────────

test_that("unweighted C++ and R agree on the sign of the log-likelihood", {
  par     <- make_par(spOccPnts)
  r_val   <- loglik_unweighted_math(spOccPnts, samMPts, par$mu, par$L)
  cpp_val <- loglik_unweighted_cpp(
    as.matrix(spOccPnts), as.matrix(samMPts), par$mu, par$L
  )
  expect_equal(sign(r_val), sign(cpp_val))
})
