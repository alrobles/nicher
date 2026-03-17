test_that("loglik_presenceonly_math returns a finite scalar", {
  par <- get_ellip_par(spOccPnts)
  val <- loglik_presenceonly_math(spOccPnts, samMPts, par$mu, par$S)

  expect_length(val, 1)
  expect_true(is.finite(val))
})

test_that("loglik_presenceonly_math result is numeric", {
  par <- get_ellip_par(spOccPnts)
  val <- loglik_presenceonly_math(spOccPnts, samMPts, par$mu, par$S)

  expect_type(val, "double")
})

test_that("loglik_presenceonly_math returns a positive value (negative log-likelihood is non-negative)", {
  par <- get_ellip_par(spOccPnts)
  val <- loglik_presenceonly_math(spOccPnts, samMPts, par$mu, par$S)

  expect_gt(val, 0)
})
