skip_if_not_installed("terra")

env_from_m <- function(env_m, var_names, nr = 30L, nc = 30L, seed = 1L) {
  set.seed(seed)
  layers <- lapply(seq_len(ncol(env_m)), function(k) {
    terra::rast(nrows = nr, ncols = nc,
                vals = sample(env_m[, k], nr * nc, replace = TRUE))
  })
  out <- do.call(c, layers)
  names(out) <- var_names
  out
}

test_that("predict.nicher returns a one-layer SpatRaster in (0, 1]", {
  data(example_env_occ_2d, package = "nicher")
  data(example_env_m_2d,   package = "nicher")
  set.seed(42)
  fit <- optimize_niche(
    env_occ    = example_env_occ_2d,
    env_m      = example_env_m_2d,
    num_starts = 4L,
    likelihood = "weighted",
    verbose    = FALSE
  )
  expect_identical(fit$var_names, colnames(example_env_occ_2d))
  env <- env_from_m(example_env_m_2d, fit$var_names)
  s <- predict(fit, env)
  expect_s4_class(s, "SpatRaster")
  expect_equal(terra::nlyr(s), 1L)
  v <- terra::values(s, na.rm = TRUE)
  expect_true(all(v > 0 & v <= 1))
})

test_that("predict.nicher reorders env layers by name", {
  data(example_env_occ_2d, package = "nicher")
  data(example_env_m_2d,   package = "nicher")
  set.seed(7)
  fit <- optimize_niche(
    env_occ    = example_env_occ_2d,
    env_m      = example_env_m_2d,
    num_starts = 3L,
    likelihood = "presence_only",
    verbose    = FALSE
  )
  env <- env_from_m(example_env_m_2d, fit$var_names, seed = 11L)
  env_rev <- env[[rev(fit$var_names)]]
  s     <- predict(fit, env)
  s_rev <- predict(fit, env_rev)
  expect_equal(terra::values(s), terra::values(s_rev))
})

test_that("predict.nicher errors when env is missing required layer names", {
  data(example_env_occ_2d, package = "nicher")
  data(example_env_m_2d,   package = "nicher")
  set.seed(8)
  fit <- optimize_niche(
    env_occ    = example_env_occ_2d,
    env_m      = example_env_m_2d,
    num_starts = 3L,
    likelihood = "weighted",
    verbose    = FALSE
  )
  env_bad <- env_from_m(example_env_m_2d, fit$var_names, seed = 12L)
  names(env_bad) <- c("foo", "bar")
  expect_error(predict(fit, env_bad), "missing required layer")
})

test_that("predict.nicher uses the eta stored on the fit (regression)", {
  data(example_env_occ_2d, package = "nicher")
  data(example_env_m_2d,   package = "nicher")
  # Fit with a non-default eta. cvine_cholesky maps `v` to a correlation
  # Cholesky factor through a Beta-quantile mapping parameterised by eta,
  # so the same theta with different eta yields a different Sigma. If
  # predict.nicher silently used eta = 1, the suitability map would not
  # match the model that was fitted.
  set.seed(13)
  fit <- optimize_niche(
    env_occ    = example_env_occ_2d,
    env_m      = example_env_m_2d,
    num_starts = 3L,
    likelihood = "weighted",
    eta        = 2.5,
    verbose    = FALSE
  )
  expect_equal(fit$eta, 2.5)

  # Reconstruct the expected Sigma using the same eta the optimizer saw.
  theta  <- fit$best$theta
  p      <- length(fit$var_names)
  v      <- if (p > 1L) theta[(2L * p + 1L):length(theta)] else numeric(0)
  sigma  <- exp(theta[(p + 1L):(2L * p)])
  L_corr_correct <- cvine_cholesky(v, d = p, eta = 2.5)
  Sigma_correct  <- tcrossprod(diag(sigma, p) %*% L_corr_correct)

  L_corr_wrong   <- cvine_cholesky(v, d = p, eta = 1.0)
  Sigma_wrong    <- tcrossprod(diag(sigma, p) %*% L_corr_wrong)

  # Sanity: the two reconstructions actually differ — otherwise this
  # test would not be exercising the fix.
  expect_false(isTRUE(all.equal(Sigma_correct, Sigma_wrong)))

  env <- env_from_m(example_env_m_2d, fit$var_names, seed = 21L)
  s_predicted <- predict(fit, env)

  # Build the reference suitability map directly with the correct Sigma.
  ref_mat <- as.matrix(terra::values(env))
  L_inv <- backsolve(t(chol(Sigma_correct)), diag(p), upper.tri = FALSE)
  ref <- niche_suitability_cpp(
    as.numeric(ref_mat), c(nrow(ref_mat), p),
    theta[seq_len(p)], L_inv, FALSE, 0L
  )
  expect_equal(as.numeric(terra::values(s_predicted)), ref,
               tolerance = 1e-12)
})

test_that("suitability evaluated at the fitted mu equals 1", {
  data(example_env_occ_2d, package = "nicher")
  data(example_env_m_2d,   package = "nicher")
  set.seed(9)
  fit <- optimize_niche(
    env_occ    = example_env_occ_2d,
    env_m      = example_env_m_2d,
    num_starts = 3L,
    likelihood = "weighted",
    verbose    = FALSE
  )
  theta  <- fit$best$theta
  p      <- length(fit$var_names)
  mu_hat <- theta[seq_len(p)]
  sigma  <- exp(theta[(p + 1L):(2L * p)])
  v      <- if (p > 1L) theta[(2L * p + 1L):length(theta)] else numeric(0)
  L_corr <- cvine_cholesky(v, d = p, eta = 1)
  Sigma  <- tcrossprod(diag(sigma, p) %*% L_corr)
  L_inv  <- backsolve(t(chol(Sigma)), diag(p), upper.tri = FALSE)
  s_at_mu <- niche_suitability_cpp(mu_hat, c(1L, p), mu_hat, L_inv, FALSE, 0L)
  expect_equal(s_at_mu, 1.0, tolerance = 1e-12)
})
