test_that("loglik_presence_math returns a finite scalar", {
  par <- get_ellip_par(spOccPnts)
  val <- loglik_presence_math(spOccPnts, samMPts, par$mu, par$S)

  expect_length(val, 1)
  expect_true(is.finite(val))
})

test_that("loglik_presence_math result is numeric", {
  par <- get_ellip_par(spOccPnts)
  val <- loglik_presence_math(spOccPnts, samMPts, par$mu, par$S)

  expect_type(val, "double")
})

test_that("negloglike_multivariable delegates to loglik_presence_math", {
  par  <- get_ellip_par(spOccPnts)
  new  <- loglik_presence_math(spOccPnts, samMPts, par$mu, par$S)
  old  <- negloglike_multivariable(par$mu, par$S, spOccPnts, samMPts)

  expect_equal(old, new, tolerance = 1e-12)
})

test_that("loglik_presence_math returns a positive value (negative log-likelihood is non-negative)", {
  par <- get_ellip_par(spOccPnts)
  val <- loglik_presence_math(spOccPnts, samMPts, par$mu, par$S)

  expect_gt(val, 0)
})
