# tests/testthat/test-cran-readiness.R
#
# Additional tests for CRAN readiness covering:
#   loglik_niche, cvine_cholesky, kde_gaussian, start_theta / start_theta_multiple,
#   get_ellipsoid_pars, NaN/Inf penalty guards, 2D paths, optimize_niche error handling.

# ---------------------------------------------------------------------------
# 6.1 Tests for loglik_niche() (pure R baseline)
# ---------------------------------------------------------------------------

test_that("loglik_niche returns a finite scalar on 2D data", {
  val <- loglik_niche(
    mu      = example_mu_vec,
    s_mat   = example_s_mat,
    env_occ = example_env_occ_2d,
    env_m   = example_env_m_2d
  )
  expect_true(is.numeric(val) && length(val) == 1L)
  expect_true(is.finite(val))
})

test_that("loglik_niche neg=FALSE returns the negative of neg=TRUE", {
  val_neg <- loglik_niche(
    mu = example_mu_vec, s_mat = example_s_mat,
    env_occ = example_env_occ_2d, env_m = example_env_m_2d,
    neg = TRUE
  )
  val_pos <- loglik_niche(
    mu = example_mu_vec, s_mat = example_s_mat,
    env_occ = example_env_occ_2d, env_m = example_env_m_2d,
    neg = FALSE
  )
  expect_equal(val_neg, -val_pos)
})

# ---------------------------------------------------------------------------
# 6.2 Tests for cvine_cholesky()
# ---------------------------------------------------------------------------

test_that("cvine_cholesky d=1 returns a 1x1 identity", {
  L <- cvine_cholesky(numeric(0), d = 1, eta = 1)
  expect_equal(dim(L), c(1L, 1L))
  expect_equal(L[1, 1], 1)
})

test_that("cvine_cholesky d=2 gives a valid correlation Cholesky factor", {
  L <- cvine_cholesky(0.5, d = 2, eta = 1)
  R <- tcrossprod(L)
  expect_equal(diag(R), c(1, 1), tolerance = 1e-10)
  expect_true(all(eigen(R, only.values = TRUE)$values > 0))
})

test_that("cvine_cholesky d=3 gives a valid correlation Cholesky factor", {
  L <- cvine_cholesky(c(0.1, -0.2, 0.8), d = 3, eta = 1)
  R <- tcrossprod(L)
  expect_equal(diag(R), rep(1, 3), tolerance = 1e-10)
  expect_true(all(eigen(R, only.values = TRUE)$values > 0))
})

test_that("cvine_cholesky raises error when v has wrong length", {
  expect_error(cvine_cholesky(c(0.1, 0.2), d = 2, eta = 1), regexp = "length")
})

# ---------------------------------------------------------------------------
# 6.3 Tests for kde_gaussian()
# ---------------------------------------------------------------------------

test_that("kde_gaussian 2D path returns correct-length positive vector", {
  x    <- as.matrix(example_env_occ_2d)
  data <- as.matrix(example_env_m_2d)
  dens <- kde_gaussian(x, data)
  expect_length(dens, nrow(x))
  expect_true(all(dens > 0))
})

test_that("kde_gaussian p>2 path returns correct-length positive vector", {
  x    <- as.matrix(example_env_occ_3d)
  data <- as.matrix(example_env_m_3d)
  dens <- kde_gaussian(x, data)
  expect_length(dens, nrow(x))
  expect_true(all(dens > 0))
})

test_that("kde_gaussian raises error when column counts differ", {
  x    <- as.matrix(example_env_occ_2d)
  data <- as.matrix(example_env_m_3d)
  expect_error(kde_gaussian(x, data), regexp = "same number of columns")
})

# ---------------------------------------------------------------------------
# 6.4 Tests for start_theta() and start_theta_multiple()
# ---------------------------------------------------------------------------

test_that("start_theta output has correct length for 2D data", {
  p  <- ncol(example_env_occ_2d)
  th <- start_theta(example_env_occ_2d)
  expect_length(th, 2L * p + p * (p - 1L) / 2L)
  expect_true(all(is.finite(th)))
})

test_that("start_theta output has correct length for 3D data", {
  p  <- ncol(example_env_occ_3d)
  th <- start_theta(example_env_occ_3d)
  expect_length(th, 2L * p + p * (p - 1L) / 2L)
  expect_true(all(is.finite(th)))
})

test_that("start_theta_multiple uniform method produces finite starts", {
  set.seed(7)
  starts <- start_theta_multiple(
    env_data   = example_env_occ_2d,
    num_starts = 8L,
    method     = "uniform"
  )
  expect_equal(nrow(starts), 8L)
  expect_true(all(is.finite(as.matrix(starts))))
})

test_that("start_theta_multiple sobol method requires pomp", {
  skip_if_not_installed("pomp")
  set.seed(7)
  starts <- start_theta_multiple(
    env_data   = example_env_occ_2d,
    num_starts = 8L,
    method     = "sobol"
  )
  expect_equal(nrow(starts), 8L)
  expect_true(all(is.finite(as.matrix(starts))))
})

# ---------------------------------------------------------------------------
# 6.5 Tests for get_ellipsoid_pars()
# ---------------------------------------------------------------------------

test_that("get_ellipsoid_pars returns mu and s_mat of correct dimensions", {
  res <- get_ellipsoid_pars(example_env_occ_2d)
  p   <- ncol(example_env_occ_2d)
  expect_type(res, "list")
  expect_named(res, c("mu", "s_mat"))
  expect_length(res$mu, p)
  expect_equal(dim(res$s_mat), c(p, p))
})

# ---------------------------------------------------------------------------
# 6.6 NaN/Inf penalty guard test
# ---------------------------------------------------------------------------

test_that("loglik_niche_math_cpp returns finite value for extreme parameters", {
  p      <- ncol(example_env_occ_2d)
  mu_ext <- rep(1e15, p)
  ls_ext <- rep(30, p)   # sigma = exp(30), very large
  v_ext  <- rep(0, p * (p - 1L) / 2L)
  theta_ext <- c(mu_ext, ls_ext, v_ext)

  val <- loglik_niche_math_cpp(
    theta   = theta_ext,
    env_occ = example_env_occ_2d,
    env_m   = example_env_m_2d,
    eta     = 1,
    neg     = TRUE
  )
  expect_true(is.finite(val),
    label = "loglik_niche_math_cpp returns finite value (penalty) for extreme params"
  )
})

# ---------------------------------------------------------------------------
# 6.7 2D path in wrapper tests
# ---------------------------------------------------------------------------

test_that("niche_unweighted works on 2D data", {
  occ2   <- as.matrix(example_env_occ_2d)
  M2     <- as.matrix(example_env_m_2d)
  theta2 <- start_theta(example_env_occ_2d)

  res <- niche_unweighted(occ = occ2, M = M2, start = theta2)
  expect_type(res, "list")
  expect_true(is.finite(res$value), label = "niche_unweighted 2D value is finite")
})

test_that("niche_presence_only works on 2D data", {
  occ2   <- as.matrix(example_env_occ_2d)
  theta2 <- start_theta(example_env_occ_2d)

  res <- niche_presence_only(occ = occ2, start = theta2)
  expect_type(res, "list")
  expect_true(is.finite(res$value), label = "niche_presence_only 2D value is finite")
})

# ---------------------------------------------------------------------------
# 6.8 optimize_niche() error handling
# ---------------------------------------------------------------------------

test_that("optimize_niche errors when env_m is missing for unweighted model", {
  expect_error(
    optimize_niche(
      env_occ    = example_env_occ_2d,
      env_m      = NULL,
      num_starts = 2L,
      breadth    = 0.1,
      likelihood = "unweighted"
    ),
    regexp = "env_m must be provided"
  )
})

test_that("optimize_niche errors when env_occ and env_m have different column names", {
  occ_bad <- example_env_occ_2d
  colnames(occ_bad) <- c("X1", "X2")

  expect_error(
    optimize_niche(
      env_occ    = occ_bad,
      env_m      = example_env_m_2d,
      num_starts = 2L,
      breadth    = 0.1,
      likelihood = "unweighted"
    ),
    regexp = "same variables"
  )
})
