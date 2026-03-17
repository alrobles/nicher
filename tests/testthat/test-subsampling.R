# Tests: subsampling vs full evaluation
#
# Subsampling a subset of the background points and computing the log-likelihood
# should still yield a finite scalar.  When we subsample the full dataset (100 %),
# the result must equal the full evaluation exactly.

# ── helpers ─────────────────────────────────────────────────────────────────

subsample_bg <- function(sam2, fraction, seed = 1L) {
  set.seed(seed)
  n    <- nrow(sam2)
  idx  <- sample.int(n, size = max(1L, floor(fraction * n)), replace = FALSE)
  as.matrix(sam2)[idx, , drop = FALSE]
}

# ── 100 % subsample equals full evaluation ───────────────────────────────────

test_that("unweighted_math: 100-pct subsample equals full evaluation", {
  par <- get_ellip_par(spOccPnts)
  L   <- t(chol(par$S))

  full_bg  <- as.matrix(samMPts)
  sub_bg   <- subsample_bg(samMPts, fraction = 1.0, seed = 42L)

  # Both use the same rows (just potentially reordered); log-likelihood is
  # order-invariant, so the values must be identical.
  full_val <- loglik_unweighted_math(spOccPnts, full_bg, par$mu, L)
  sub_val  <- loglik_unweighted_math(spOccPnts, sub_bg,  par$mu, L)

  expect_equal(full_val, sub_val, tolerance = 1e-10)
})

test_that("unweighted_cpp: 100-pct subsample equals full evaluation", {
  par <- get_ellip_par(spOccPnts)
  L   <- t(chol(par$S))

  full_bg <- as.matrix(samMPts)
  sub_bg  <- subsample_bg(samMPts, fraction = 1.0, seed = 42L)

  full_val <- loglik_unweighted_cpp(as.matrix(spOccPnts), full_bg, par$mu, L)
  sub_val  <- loglik_unweighted_cpp(as.matrix(spOccPnts), sub_bg,  par$mu, L)

  expect_equal(full_val, sub_val, tolerance = 1e-10)
})

test_that("presenceonly_math: 100-pct subsample equals full evaluation", {
  par     <- get_ellip_par(spOccPnts)
  full_bg <- as.matrix(samMPts)
  sub_bg  <- subsample_bg(samMPts, fraction = 1.0, seed = 99L)

  full_val <- loglik_presenceonly_math(spOccPnts, full_bg, par$mu, par$S)
  sub_val  <- loglik_presenceonly_math(spOccPnts, sub_bg,  par$mu, par$S)

  expect_equal(full_val, sub_val, tolerance = 1e-10)
})

# ── partial subsamples return a finite scalar ─────────────────────────────────

test_that("unweighted_math: 50-pct subsample returns finite scalar", {
  par     <- get_ellip_par(spOccPnts)
  L       <- t(chol(par$S))
  sub_bg  <- subsample_bg(samMPts, fraction = 0.5, seed = 1L)

  val <- loglik_unweighted_math(spOccPnts, sub_bg, par$mu, L)

  expect_length(val, 1L)
  expect_true(is.finite(val))
})

test_that("unweighted_cpp: 50-pct subsample returns finite scalar", {
  par     <- get_ellip_par(spOccPnts)
  L       <- t(chol(par$S))
  sub_bg  <- subsample_bg(samMPts, fraction = 0.5, seed = 2L)

  val <- loglik_unweighted_cpp(as.matrix(spOccPnts), sub_bg, par$mu, L)

  expect_length(val, 1L)
  expect_true(is.finite(val))
})

test_that("presenceonly_math: 50-pct subsample returns finite scalar", {
  par    <- get_ellip_par(spOccPnts)
  sub_bg <- subsample_bg(samMPts, fraction = 0.5, seed = 3L)

  val <- loglik_presenceonly_math(spOccPnts, sub_bg, par$mu, par$S)

  expect_length(val, 1L)
  expect_true(is.finite(val))
})

test_that("presenceonly_cpp: 50-pct subsample returns finite scalar", {
  par    <- get_ellip_par(spOccPnts)
  sub_bg <- subsample_bg(samMPts, fraction = 0.5, seed = 4L)

  val <- loglik_presenceonly_cpp(
    as.matrix(spOccPnts), sub_bg, par$mu, par$S
  )

  expect_length(val, 1L)
  expect_true(is.finite(val))
})

# ── R and C++ agree on subsampled data ───────────────────────────────────────

test_that("unweighted R and C++ agree on 50-pct subsampled background", {
  par    <- get_ellip_par(spOccPnts)
  L      <- t(chol(par$S))
  sub_bg <- subsample_bg(samMPts, fraction = 0.5, seed = 7L)

  r_val   <- loglik_unweighted_math(spOccPnts, sub_bg, par$mu, L)
  cpp_val <- loglik_unweighted_cpp(as.matrix(spOccPnts), sub_bg, par$mu, L)

  expect_equal(r_val, cpp_val, tolerance = 1e-8)
})

test_that("presenceonly R and C++ agree on 50-pct subsampled background", {
  par    <- get_ellip_par(spOccPnts)
  sub_bg <- subsample_bg(samMPts, fraction = 0.5, seed = 8L)

  r_val   <- loglik_presenceonly_math(spOccPnts, sub_bg, par$mu, par$S)
  cpp_val <- loglik_presenceonly_cpp(as.matrix(spOccPnts), sub_bg, par$mu, par$S)

  expect_equal(r_val, cpp_val, tolerance = 1e-8)
})

# ── larger vs smaller subsamples differ ──────────────────────────────────────

test_that("different subsample sizes give different unweighted log-likelihoods", {
  par     <- get_ellip_par(spOccPnts)
  L       <- t(chol(par$S))
  sub25   <- subsample_bg(samMPts, fraction = 0.25, seed = 11L)
  sub75   <- subsample_bg(samMPts, fraction = 0.75, seed = 11L)

  val25 <- loglik_unweighted_math(spOccPnts, sub25, par$mu, L)
  val75 <- loglik_unweighted_math(spOccPnts, sub75, par$mu, L)

  expect_false(isTRUE(all.equal(val25, val75, tolerance = 1e-6)))
})

# ── single-point background ───────────────────────────────────────────────────

test_that("presenceonly_math handles single background point", {
  par    <- get_ellip_par(spOccPnts)
  one_bg <- as.matrix(samMPts)[1L, , drop = FALSE]

  val <- loglik_presenceonly_math(spOccPnts, one_bg, par$mu, par$S)

  expect_length(val, 1L)
  expect_true(is.finite(val))
})

test_that("presenceonly_cpp handles single background point", {
  par    <- get_ellip_par(spOccPnts)
  one_bg <- as.matrix(samMPts)[1L, , drop = FALSE]

  val <- loglik_presenceonly_cpp(
    as.matrix(spOccPnts), one_bg, par$mu, par$S
  )

  expect_length(val, 1L)
  expect_true(is.finite(val))
})
