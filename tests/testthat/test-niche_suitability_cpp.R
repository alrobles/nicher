test_that("centroid suitability is exactly 1", {
  set.seed(0)
  for (p in c(1L, 2L, 4L)) {
    A <- matrix(stats::rnorm(p * p), p, p)
    Sigma <- crossprod(A) + diag(p)
    mu <- stats::runif(p)
    L <- t(chol(Sigma))
    L_inv <- backsolve(L, diag(p), upper.tri = FALSE)
    out <- niche_suitability_cpp(
      env_dat_vec  = mu,
      env_dat_dims = c(1L, p),
      mu           = mu,
      L_inv        = L_inv,
      return_log   = FALSE,
      num_threads  = 0L
    )
    expect_equal(out, 1.0, tolerance = 1e-12)
  }
})

test_that("vectorized suitability matches stats::mahalanobis to machine precision", {
  set.seed(1)
  p <- 4L
  n <- 200L
  A <- matrix(stats::rnorm(p * p), p, p)
  Sigma <- crossprod(A) + diag(p)
  mu <- stats::rnorm(p)
  X <- matrix(stats::rnorm(n * p), n, p)
  L <- t(chol(Sigma))
  L_inv <- backsolve(L, diag(p), upper.tri = FALSE)

  ref <- exp(-0.5 * stats::mahalanobis(X, center = mu, cov = Sigma))
  got <- niche_suitability_cpp(
    env_dat_vec  = as.numeric(X),
    env_dat_dims = c(n, p),
    mu           = mu,
    L_inv        = L_inv,
    return_log   = FALSE,
    num_threads  = 0L
  )
  expect_equal(got, ref, tolerance = 1e-12)
  expect_true(all(got > 0 & got <= 1))
})

test_that("return_log = TRUE round-trips with exp() to machine precision", {
  set.seed(2)
  p <- 3L; n <- 100L
  A <- matrix(stats::rnorm(p * p), p, p); Sigma <- crossprod(A) + diag(p)
  mu <- stats::rnorm(p); X <- matrix(stats::rnorm(n * p), n, p)
  L <- t(chol(Sigma)); L_inv <- backsolve(L, diag(p), upper.tri = FALSE)

  s <- niche_suitability_cpp(as.numeric(X), c(n, p), mu, L_inv,
                             return_log = FALSE, num_threads = 0L)
  l <- niche_suitability_cpp(as.numeric(X), c(n, p), mu, L_inv,
                             return_log = TRUE,  num_threads = 0L)
  expect_equal(exp(l), s, tolerance = 1e-14)
  expect_true(all(l <= 0))
})

test_that("NA inputs produce NA outputs, not aborts", {
  set.seed(3)
  p <- 2L; n <- 5L
  X <- matrix(stats::rnorm(n * p), n, p)
  X[2, 1] <- NA_real_
  X[5, 2] <- NaN
  mu <- c(0, 0); Sigma <- diag(p)
  L_inv <- backsolve(t(chol(Sigma)), diag(p), upper.tri = FALSE)
  out <- niche_suitability_cpp(as.numeric(X), c(n, p), mu, L_inv, FALSE, 0L)
  expect_true(is.na(out[2]))
  expect_true(is.na(out[5]))
  expect_true(all(is.finite(out[c(1, 3, 4)])))
})

test_that("input validation rejects malformed args", {
  expect_error(
    niche_suitability_cpp(numeric(4), c(2L), c(0, 0),
                          diag(2L), FALSE, 0L),
    "env_dat_dims"
  )
  expect_error(
    niche_suitability_cpp(numeric(3), c(2L, 2L), c(0, 0),
                          diag(2L), FALSE, 0L),
    "n_loc \\* p"
  )
  expect_error(
    niche_suitability_cpp(numeric(4), c(2L, 2L), c(0, 0, 0),
                          diag(2L), FALSE, 0L),
    "length p"
  )
  expect_error(
    niche_suitability_cpp(numeric(4), c(2L, 2L), c(0, 0),
                          diag(3L), FALSE, 0L),
    "p x p"
  )
})

test_that("results are deterministic across thread counts", {
  set.seed(4)
  p <- 3L; n <- 5000L
  A <- matrix(stats::rnorm(p * p), p, p); Sigma <- crossprod(A) + diag(p)
  mu <- stats::rnorm(p); X <- matrix(stats::rnorm(n * p), n, p)
  L_inv <- backsolve(t(chol(Sigma)), diag(p), upper.tri = FALSE)
  out1 <- niche_suitability_cpp(as.numeric(X), c(n, p), mu, L_inv, FALSE, 1L)
  out4 <- niche_suitability_cpp(as.numeric(X), c(n, p), mu, L_inv, FALSE, 4L)
  expect_identical(out1, out4)
})
