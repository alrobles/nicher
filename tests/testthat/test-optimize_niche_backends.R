test_that("optimize_niche backend='cpp' and 'r' agree on the same Sobol starts", {
  occ <- example_env_occ_2d
  M   <- example_env_m_2d

  fit_cpp <- suppressWarnings(optimize_niche(
    env_occ = occ, env_m = M,
    num_starts = 5L, breadth = 0.1,
    likelihood = "weighted",
    backend    = "cpp",
    seed       = 1L,
    control    = list(maxeval = 500L)
  ))

  fit_r <- suppressWarnings(optimize_niche(
    env_occ = occ, env_m = M,
    num_starts = 5L, breadth = 0.1,
    likelihood = "weighted",
    backend    = "r",
    seed       = 1L,
    control    = list(maxeval = 500L)
  ))

  expect_s3_class(fit_cpp, "nicher")
  expect_s3_class(fit_r,   "nicher")

  # The two backends should reach a similar best log-likelihood. We allow a
  # generous tolerance because the C++ backend uses ucminfcpp's BFGS while
  # the R backend uses ucminf's BFGS, and FD step choices differ slightly.
  expect_equal(fit_cpp$best$loglik, fit_r$best$loglik, tolerance = 1e-2)
})

test_that("backend='r' emits a lifecycle deprecation warning", {
  occ <- example_env_occ_2d
  M   <- example_env_m_2d
  expect_warning(
    optimize_niche(
      env_occ = occ, env_m = M,
      num_starts = 2L, likelihood = "weighted",
      backend = "r", seed = 1L,
      control = list(maxeval = 100L)
    ),
    regexp = "deprecated|backend"
  )
})
