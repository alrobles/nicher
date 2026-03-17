test_that("loglik_weighted_math returns a finite scalar with default weights", {
  par <- get_ellip_par(spOccPnts)
  val <- loglik_weighted_math(spOccPnts, samMPts, par$mu, par$S)

  expect_length(val, 1)
  expect_true(is.finite(val))
})

test_that("loglik_weighted_math result is numeric", {
  par <- get_ellip_par(spOccPnts)
  val <- loglik_weighted_math(spOccPnts, samMPts, par$mu, par$S)

  expect_type(val, "double")
})

test_that("loglik_weighted_math with uniform weights equals loglik_presenceonly_math", {
  par <- get_ellip_par(spOccPnts)
  val_uniform <- loglik_weighted_math(
    spOccPnts, samMPts, par$mu, par$S,
    weights = rep(1, nrow(samMPts))
  )
  val_presence <- loglik_presenceonly_math(spOccPnts, samMPts, par$mu, par$S)

  expect_equal(val_uniform, val_presence, tolerance = 1e-10)
})

test_that("loglik_weighted_math with NULL weights equals loglik_presenceonly_math", {
  par <- get_ellip_par(spOccPnts)
  val_null     <- loglik_weighted_math(spOccPnts, samMPts, par$mu, par$S)
  val_presence <- loglik_presenceonly_math(spOccPnts, samMPts, par$mu, par$S)

  expect_equal(val_null, val_presence, tolerance = 1e-10)
})

test_that("loglik_weighted_math weights are normalized invariant to scale", {
  par <- get_ellip_par(spOccPnts)
  w   <- seq_len(nrow(samMPts))
  val_w    <- loglik_weighted_math(spOccPnts, samMPts, par$mu, par$S, weights = w)
  val_2w   <- loglik_weighted_math(spOccPnts, samMPts, par$mu, par$S, weights = 2 * w)

  expect_equal(val_w, val_2w, tolerance = 1e-12)
})

test_that("loglik_weighted_math errors on wrong-length weights", {
  par <- get_ellip_par(spOccPnts)
  expect_error(
    loglik_weighted_math(spOccPnts, samMPts, par$mu, par$S, weights = c(1, 2)),
    "'weights' must have the same length"
  )
})

test_that("loglik_weighted_math errors on negative weights", {
  par <- get_ellip_par(spOccPnts)
  w   <- rep(1, nrow(samMPts))
  w[1] <- -1
  expect_error(
    loglik_weighted_math(spOccPnts, samMPts, par$mu, par$S, weights = w),
    "'weights' must be non-negative"
  )
})

test_that("loglik_weighted_math errors on all-zero weights", {
  par <- get_ellip_par(spOccPnts)
  expect_error(
    loglik_weighted_math(spOccPnts, samMPts, par$mu, par$S,
                         weights = rep(0, nrow(samMPts))),
    "'weights' must not all be zero"
  )
})
