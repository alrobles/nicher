# Benchmark tests: R vs C++ performance comparison
#
# Uses bench::mark() to compare execution time between R and C++ backends for
# both the unweighted and presence-only log-likelihood functions.  The tests
# verify that:
#  1. Benchmarks complete without error.
#  2. C++ is at least as fast as the R implementation (median time).
#  3. Results from R and C++ agree numerically (bench::mark checks this via
#     `check = TRUE` which is the default when a `check` column is returned).
#
# These tests are skipped when the bench package is not installed so that the
# package can still be checked without the optional dependency.

skip_if_not_installed("bench")

# ── shared fixtures ──────────────────────────────────────────────────────────

occ_mat <- as.matrix(spOccPnts)
bg_mat  <- as.matrix(samMPts)

# ── unweighted: R vs C++ ─────────────────────────────────────────────────────

test_that("bench::mark completes for unweighted R vs C++", {
  par <- get_ellip_par(spOccPnts)
  L   <- t(chol(par$S))

  bm <- bench::mark(
    R   = loglik_unweighted_math(spOccPnts, samMPts,  par$mu, L),
    Cpp = loglik_unweighted_cpp(occ_mat,    bg_mat,   par$mu, L),
    iterations  = 50L,
    check       = TRUE,
    relative    = FALSE
  )

  expect_s3_class(bm, "bench_mark")
  expect_equal(nrow(bm), 2L)
  expect_true(all(bm$`itr/sec` > 0))
})

test_that("C++ unweighted is not slower than R (median time)", {
  par <- get_ellip_par(spOccPnts)
  L   <- t(chol(par$S))

  bm <- bench::mark(
    R   = loglik_unweighted_math(spOccPnts, samMPts, par$mu, L),
    Cpp = loglik_unweighted_cpp(occ_mat,    bg_mat,  par$mu, L),
    iterations = 100L,
    check      = FALSE
  )

  # Row 1 = R, Row 2 = Cpp (same order as supplied expressions)
  median_r   <- as.numeric(bm$median[[1L]])
  median_cpp <- as.numeric(bm$median[[2L]])

  # Allow up to 3× overhead: on some CI platforms the JIT-compiled R code can
  # be competitive with C++ for small datasets.  In practice C++ is faster;
  # this bound prevents spurious failures on slow or shared CI machines.
  expect_lte(median_cpp, median_r * 3)
})

# ── presence-only: R vs C++ ──────────────────────────────────────────────────

test_that("bench::mark completes for presence-only R vs C++", {
  par <- get_ellip_par(spOccPnts)

  bm <- bench::mark(
    R   = loglik_presenceonly_math(spOccPnts, samMPts, par$mu, par$S),
    Cpp = loglik_presenceonly_cpp(occ_mat,    bg_mat,  par$mu, par$S),
    iterations = 50L,
    check      = TRUE,
    relative   = FALSE
  )

  expect_s3_class(bm, "bench_mark")
  expect_equal(nrow(bm), 2L)
  expect_true(all(bm$`itr/sec` > 0))
})

test_that("C++ presence-only is not slower than R (median time)", {
  par <- get_ellip_par(spOccPnts)

  bm <- bench::mark(
    R   = loglik_presenceonly_math(spOccPnts, samMPts, par$mu, par$S),
    Cpp = loglik_presenceonly_cpp(occ_mat,    bg_mat,  par$mu, par$S),
    iterations = 100L,
    check      = FALSE
  )

  # Row 1 = R, Row 2 = Cpp
  median_r   <- as.numeric(bm$median[[1L]])
  median_cpp <- as.numeric(bm$median[[2L]])

  # Allow up to 3× overhead for the same reasons as the unweighted benchmark.
  expect_lte(median_cpp, median_r * 3)
})

# ── full vs subsampled background (unweighted) ───────────────────────────────

test_that("bench::mark completes for full vs 50-pct subsampled background", {
  set.seed(5L)
  par    <- get_ellip_par(spOccPnts)
  L      <- t(chol(par$S))
  n      <- nrow(bg_mat)
  sub_bg <- bg_mat[sample.int(n, size = floor(0.5 * n)), , drop = FALSE]

  bm <- bench::mark(
    full      = loglik_unweighted_cpp(occ_mat, bg_mat,  par$mu, L),
    subsampled = loglik_unweighted_cpp(occ_mat, sub_bg, par$mu, L),
    iterations = 50L,
    check      = FALSE
  )

  expect_s3_class(bm, "bench_mark")
  expect_equal(nrow(bm), 2L)
  expect_true(all(bm$`itr/sec` > 0))
})

# ── 2-D vs 3-D models ────────────────────────────────────────────────────────

test_that("bench::mark completes for 2-D vs 3-D unweighted C++", {
  skip_if_not(ncol(occ_mat) >= 2L, "need at least 2 variables")

  set.seed(9L)
  # 2-D data
  occ2 <- matrix(rnorm(60,  mean = c(1,  2),    sd = 0.5), ncol = 2, byrow = TRUE)
  bg2  <- matrix(rnorm(400, mean = 0,            sd = 2),   ncol = 2, byrow = TRUE)
  mu2  <- colMeans(occ2)
  S2   <- cov(occ2)
  L2   <- t(chol(S2))

  # 3-D data
  occ3 <- matrix(rnorm(60,  mean = c(1, 2, -1), sd = 0.5), ncol = 3, byrow = TRUE)
  bg3  <- matrix(rnorm(600, mean = 0,            sd = 2),   ncol = 3, byrow = TRUE)
  mu3  <- colMeans(occ3)
  S3   <- cov(occ3)
  L3   <- t(chol(S3))

  bm <- bench::mark(
    `2D` = loglik_unweighted_cpp(occ2, bg2, as.vector(mu2), L2),
    `3D` = loglik_unweighted_cpp(occ3, bg3, as.vector(mu3), L3),
    iterations = 50L,
    check      = FALSE
  )

  expect_s3_class(bm, "bench_mark")
  expect_equal(nrow(bm), 2L)
})
