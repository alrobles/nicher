# Stage 4 integration tests: verify that optimize_niche() converges for
# the NCST (Non-Central Skew t) likelihood families.

test_that("optimize_niche(likelihood = 'ncst') converges on 2D data", {
  set.seed(42)
  fit <- optimize_niche(
    env_occ    = example_env_occ_2d,
    num_starts = 2L,
    likelihood = "ncst"
  )

  expect_s3_class(fit, "nicher")
  expect_true(is.finite(fit$best$loglik))
  expect_equal(fit$likelihood, "ncst")

  # theta should have correct length: 3*p + n_v + 1
  p <- ncol(example_env_occ_2d)
  n_v <- p * (p - 1) / 2
  expect_length(fit$best$theta, 3 * p + n_v + 1)

  # predict uses (mu, Sigma) Gaussian projection
  pars <- nicher:::.recover_mu_sigma(fit)
  expect_equal(length(pars$mu), p)
  expect_equal(nrow(pars$Sigma), p)
})

test_that("optimize_niche(likelihood = 'ncst_weighted') converges on 2D data", {
  set.seed(42)
  # Use tiny synthetic data: NCST weighted has 32-node quadrature
  # applied to every M row, so large M is prohibitively slow in CI.
  p <- 2L
  small_occ <- example_env_occ_2d[sample.int(nrow(example_env_occ_2d), 30L), ]
  small_m   <- example_env_m_2d[sample.int(nrow(example_env_m_2d), 100L), ]
  fit <- suppressWarnings(optimize_niche(
    env_occ    = small_occ,
    env_m      = small_m,
    num_starts = 2L,
    likelihood = "ncst_weighted",
    warm_start = FALSE,
    prior_mu_lambda = 0,
    prior_log_sigma_lambda = 0,
    control = list(maxeval = 200L)
  ))

  expect_s3_class(fit, "nicher")
  expect_true(is.finite(fit$best$loglik))
  expect_equal(fit$likelihood, "ncst_weighted")

  n_v <- p * (p - 1) / 2
  expect_length(fit$best$theta, 3 * p + n_v + 1)
})
