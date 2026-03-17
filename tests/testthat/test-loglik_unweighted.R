test_that("loglik_unweighted_math returns a finite scalar", {
  par <- get_ellip_par(spOccPnts)
  L   <- t(chol(par$S))
  val <- loglik_unweighted_math(spOccPnts, samMPts, par$mu, L)

  expect_length(val, 1)
  expect_true(is.finite(val))
})

test_that("loglik_unweighted_math result is numeric", {
  par <- get_ellip_par(spOccPnts)
  L   <- t(chol(par$S))
  val <- loglik_unweighted_math(spOccPnts, samMPts, par$mu, L)

  expect_type(val, "double")
})

test_that("loglik_unweighted_math accepts data frame inputs", {
  par <- get_ellip_par(spOccPnts)
  L   <- t(chol(par$S))
  val_df  <- loglik_unweighted_math(spOccPnts, samMPts, par$mu, L)
  val_mat <- loglik_unweighted_math(
    as.matrix(spOccPnts), as.matrix(samMPts), par$mu, L
  )

  expect_equal(val_df, val_mat)
})

test_that("loglik_math_niche calls loglik_unweighted_math via math_to_bio_niche", {
  par <- get_ellip_par(spOccPnts)
  L   <- t(chol(par$S))
  p   <- length(par$mu)
  k   <- p * (p + 1) / 2

  # Build math-scale param_vector
  L_log       <- L
  diag(L_log) <- log(diag(L))
  param_vec   <- c(par$mu, L_log[lower.tri(L_log, diag = TRUE)])
  names(param_vec) <- c(paste0("mu", seq_len(p)), paste0("L", seq_len(k)))

  val_math  <- loglik_math_niche(param_vec, spOccPnts, samMPts, negative = FALSE)
  val_bio   <- loglik_unweighted_math(spOccPnts, samMPts, par$mu, L)

  expect_equal(val_math, val_bio, tolerance = 1e-8)
})
