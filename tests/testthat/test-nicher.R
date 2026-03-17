## Tests for the nicher() wrapper function and its S3 methods

# ---------------------------------------------------------------------------
# nicher() basic smoke tests
# ---------------------------------------------------------------------------

test_that("nicher() returns a nicher object for presence_only model", {
  fit <- nicher(spOccPnts, samMPts, model = "presence_only")
  expect_s3_class(fit, "nicher")
  expect_true(is.nicher(fit))
})

test_that("nicher() returns a nicher object for unweighted model", {
  fit <- nicher(spOccPnts, samMPts, model = "unweighted")
  expect_s3_class(fit, "nicher")
})

test_that("nicher() returns a nicher object for weighted model", {
  fit <- nicher(spOccPnts, samMPts, model = "weighted")
  expect_s3_class(fit, "nicher")
})

# ---------------------------------------------------------------------------
# nicher_object structure
# ---------------------------------------------------------------------------

test_that("nicher object has expected components", {
  fit <- nicher(spOccPnts, samMPts, model = "presence_only")

  expect_named(fit, c("loglik", "math_params", "bio_params", "optim_table", "model", "method"),
               ignore.order = TRUE)
})

test_that("nicher object stores correct method and model", {
  fit <- nicher(spOccPnts, samMPts, method = "mle", model = "presence_only")

  expect_equal(fit$method, "mle")
  expect_equal(fit$model,  "presence_only")
})

test_that("loglik is a finite scalar", {
  fit <- nicher(spOccPnts, samMPts, model = "presence_only")

  expect_length(fit$loglik, 1L)
  expect_true(is.finite(fit$loglik))
})

# ---------------------------------------------------------------------------
# math_params structure and naming
# ---------------------------------------------------------------------------

test_that("math_params has canonical names (mu1..mup, L1..Lk)", {
  fit <- nicher(spOccPnts, samMPts, model = "presence_only")
  p   <- ncol(spOccPnts)
  k   <- p * (p + 1L) / 2L

  expected_names <- c(paste0("mu", seq_len(p)), paste0("L", seq_len(k)))
  expect_equal(names(fit$math_params), expected_names)
})

test_that("math_params is a finite numeric vector", {
  fit <- nicher(spOccPnts, samMPts, model = "presence_only")

  expect_true(is.numeric(fit$math_params))
  expect_true(all(is.finite(fit$math_params)))
})

# ---------------------------------------------------------------------------
# bio_params structure
# ---------------------------------------------------------------------------

test_that("bio_params contains mu, S, sigma, R", {
  fit <- nicher(spOccPnts, samMPts, model = "presence_only")

  expect_named(fit$bio_params, c("mu", "S", "sigma", "R"), ignore.order = TRUE)
})

test_that("bio_params$mu has length p", {
  fit <- nicher(spOccPnts, samMPts, model = "presence_only")
  p   <- ncol(spOccPnts)

  expect_length(fit$bio_params$mu, p)
})

test_that("bio_params$S is a symmetric p x p matrix", {
  fit <- nicher(spOccPnts, samMPts, model = "presence_only")
  p   <- ncol(spOccPnts)
  S   <- fit$bio_params$S

  expect_equal(dim(S), c(p, p))
  expect_equal(S, t(S), tolerance = 1e-10)
})

test_that("bio_params$R is a symmetric p x p correlation matrix with unit diagonal", {
  fit <- nicher(spOccPnts, samMPts, model = "presence_only")
  p   <- ncol(spOccPnts)
  R   <- fit$bio_params$R

  expect_equal(dim(R), c(p, p))
  expect_equal(diag(R), rep(1, p), tolerance = 1e-10)
  expect_equal(R, t(R), tolerance = 1e-10)
})

test_that("bio_params$sigma equals sqrt of diagonal of S", {
  fit   <- nicher(spOccPnts, samMPts, model = "presence_only")
  sigma <- fit$bio_params$sigma
  S     <- fit$bio_params$S

  expect_equal(sigma, sqrt(diag(S)), tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# optim_table
# ---------------------------------------------------------------------------

test_that("optim_table is a data frame with at least one row", {
  fit <- nicher(spOccPnts, samMPts, model = "presence_only")

  expect_s3_class(fit$optim_table, "data.frame")
  expect_gte(nrow(fit$optim_table), 1L)
})

# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

test_that("print.nicher returns the object invisibly", {
  fit <- nicher(spOccPnts, samMPts, model = "presence_only")

  out <- capture.output(result <- print(fit))
  expect_identical(result, fit)
  expect_true(any(grepl("Nicher model", out)))
  expect_true(any(grepl("Log-likelihood", out)))
})

test_that("summary.nicher returns the object invisibly", {
  fit <- nicher(spOccPnts, samMPts, model = "presence_only")

  out <- capture.output(result <- summary(fit))
  expect_identical(result, fit)
  expect_true(any(grepl("nicher model summary", out)))
  expect_true(any(grepl("Covariance matrix", out)))
  expect_true(any(grepl("math", out, ignore.case = TRUE)))
})

# ---------------------------------------------------------------------------
# is.nicher
# ---------------------------------------------------------------------------

test_that("is.nicher returns TRUE for nicher objects and FALSE otherwise", {
  fit <- nicher(spOccPnts, samMPts, model = "presence_only")

  expect_true(is.nicher(fit))
  expect_false(is.nicher(42))
  expect_false(is.nicher(list()))
})

# ---------------------------------------------------------------------------
# nicher_object constructor
# ---------------------------------------------------------------------------

test_that("nicher_object creates a nicher S3 object", {
  par <- get_ellip_par(spOccPnts)
  L   <- t(chol(par$S))
  p   <- length(par$mu)
  k   <- p * (p + 1L) / 2L
  L_log       <- L
  diag(L_log) <- log(diag(L))
  math_params <- c(par$mu, L_log[lower.tri(L_log, diag = TRUE)])
  names(math_params) <- c(paste0("mu", seq_len(p)), paste0("L", seq_len(k)))
  S     <- L %*% t(L)
  sigma <- sqrt(diag(S))
  D_inv <- diag(1 / sigma, nrow = p)
  R     <- (D_inv %*% S %*% D_inv + t(D_inv %*% S %*% D_inv)) / 2
  bio   <- list(mu = par$mu, S = S, sigma = sigma, R = R)
  ll    <- loglik_unweighted_math(spOccPnts, samMPts, par$mu, L)

  obj <- nicher_object(ll, math_params, bio, data.frame(), "unweighted", "mle")

  expect_s3_class(obj, "nicher")
  expect_equal(obj$loglik, ll)
  expect_equal(obj$model,  "unweighted")
  expect_equal(obj$method, "mle")
})

# ---------------------------------------------------------------------------
# invalid inputs
# ---------------------------------------------------------------------------

test_that("nicher() errors on unsupported method", {
  expect_error(
    nicher(spOccPnts, samMPts, method = "bayesian"),
    "should be one of"
  )
})

test_that("nicher() errors on unknown model", {
  expect_error(
    nicher(spOccPnts, samMPts, model = "not_a_model"),
    "should be one of"
  )
})
