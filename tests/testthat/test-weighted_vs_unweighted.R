# Tests: weighted vs unweighted log-likelihood comparison
#
# The "weighted" model generalises the unweighted one: when all observation
# weights are uniform (w_i = 1/n), the weighted log-likelihood reduces to the
# unweighted formula.  These tests verify that mathematical equivalence.

# ── local weighted helper (not exported from package) ────────────────────────

#' Weighted version of the unweighted Gaussian log-likelihood.
#'
#' When `w1 = rep(1/n1, n1)` and `w2 = rep(1/n2, n2)` (uniform), the result
#' must match `loglik_unweighted_math()` exactly.
#'
#' @param sam1,sam2 Matrices of presence / background points.
#' @param mu Numeric mean vector.
#' @param L Lower-triangular Cholesky factor.
#' @param w1 Presence weights (must sum to 1); defaults to uniform.
#' @param w2 Background weights (must sum to 1); defaults to uniform.
#' @return Weighted log-likelihood (scalar).
loglik_weighted_local <- function(sam1, sam2, mu, L, w1 = NULL, w2 = NULL) {
  sam1 <- as.matrix(sam1)
  sam2 <- as.matrix(sam2)

  n1 <- nrow(sam1)
  n2 <- nrow(sam2)

  if (is.null(w1)) w1 <- rep(1 / n1, n1)
  if (is.null(w2)) w2 <- rep(1 / n2, n2)

  S      <- L %*% t(L)
  logdet <- as.numeric(determinant(S, logarithm = TRUE)$modulus)
  q1     <- stats::mahalanobis(sam1, center = mu, cov = S, inverted = FALSE)
  q2     <- stats::mahalanobis(sam2, center = mu, cov = S, inverted = FALSE)

  # Weighted analogue of the unweighted formula:
  #   logL = -0.5 * (n1 * logdet + n1 * sum(w1 * q1))
  #        + 0.5 * (n2 * logdet + n2 * sum(w2 * q2))
  logL <- -0.5 * (n1 * logdet + n1 * sum(w1 * q1)) +
           0.5  * (n2 * logdet + n2 * sum(w2 * q2))
  logL
}

# ── uniform weights ≡ unweighted ─────────────────────────────────────────────

test_that("weighted model with uniform weights equals unweighted model (package data)", {
  par <- get_ellip_par(spOccPnts)
  L   <- t(chol(par$S))
  n1  <- nrow(spOccPnts)
  n2  <- nrow(samMPts)

  unweighted <- loglik_unweighted_math(spOccPnts, samMPts, par$mu, L)
  weighted   <- loglik_weighted_local(
    spOccPnts, samMPts, par$mu, L,
    w1 = rep(1 / n1, n1),
    w2 = rep(1 / n2, n2)
  )

  expect_equal(weighted, unweighted, tolerance = 1e-10)
})

test_that("weighted model with uniform weights equals unweighted model (synthetic 2-D)", {
  set.seed(21)
  occ <- matrix(rnorm(40, mean = c(1, 2), sd = 0.5), ncol = 2, byrow = TRUE)
  bg  <- matrix(rnorm(200, mean = 0, sd = 2),          ncol = 2, byrow = TRUE)
  mu  <- colMeans(occ)
  S   <- cov(occ)
  L   <- t(chol(S))
  n1  <- nrow(occ)
  n2  <- nrow(bg)

  unweighted <- loglik_unweighted_math(occ, bg, mu, L)
  weighted   <- loglik_weighted_local(
    occ, bg, mu, L,
    w1 = rep(1 / n1, n1),
    w2 = rep(1 / n2, n2)
  )

  expect_equal(weighted, unweighted, tolerance = 1e-10)
})

test_that("weighted model with uniform weights equals unweighted model (synthetic 3-D)", {
  set.seed(77)
  occ <- matrix(rnorm(60, mean = c(0, 1, -1), sd = 0.4), ncol = 3, byrow = TRUE)
  bg  <- matrix(rnorm(300, mean = 0, sd = 2),              ncol = 3, byrow = TRUE)
  mu  <- colMeans(occ)
  S   <- cov(occ)
  L   <- t(chol(S))
  n1  <- nrow(occ)
  n2  <- nrow(bg)

  unweighted <- loglik_unweighted_math(occ, bg, mu, L)
  weighted   <- loglik_weighted_local(
    occ, bg, mu, L,
    w1 = rep(1 / n1, n1),
    w2 = rep(1 / n2, n2)
  )

  expect_equal(weighted, unweighted, tolerance = 1e-10)
})

# ── non-uniform weights alter the log-likelihood ─────────────────────────────

test_that("non-uniform weights differ from unweighted result", {
  set.seed(55)
  occ <- matrix(rnorm(40, mean = c(1, 2), sd = 0.5), ncol = 2, byrow = TRUE)
  bg  <- matrix(rnorm(200, mean = 0, sd = 2),          ncol = 2, byrow = TRUE)
  mu  <- colMeans(occ)
  S   <- cov(occ)
  L   <- t(chol(S))
  n1  <- nrow(occ)
  n2  <- nrow(bg)

  # Upweight the first presence point
  w1       <- rep(1 / n1, n1)
  w1[1]    <- w1[1] * 5
  w1       <- w1 / sum(w1)   # renormalise

  unweighted <- loglik_unweighted_math(occ, bg, mu, L)
  weighted   <- loglik_weighted_local(occ, bg, mu, L, w1 = w1,
                                      w2 = rep(1 / n2, n2))

  expect_false(isTRUE(all.equal(weighted, unweighted, tolerance = 1e-6)))
})

# ── weighted result is always a finite scalar ─────────────────────────────────

test_that("weighted local helper returns a finite scalar", {
  par <- get_ellip_par(spOccPnts)
  L   <- t(chol(par$S))
  val <- loglik_weighted_local(spOccPnts, samMPts, par$mu, L)

  expect_length(val, 1L)
  expect_true(is.finite(val))
  expect_type(val, "double")
})

# ── weighted C++ (unweighted_cpp) matches weighted R under uniform weights ────

test_that("unweighted_cpp equals weighted local with uniform weights", {
  set.seed(33)
  occ <- matrix(rnorm(50, mean = c(1, 0), sd = 0.5), ncol = 2, byrow = TRUE)
  bg  <- matrix(rnorm(250, mean = 0, sd = 2),          ncol = 2, byrow = TRUE)
  mu  <- colMeans(occ)
  S   <- cov(occ)
  L   <- t(chol(S))
  n1  <- nrow(occ)
  n2  <- nrow(bg)

  cpp_val  <- loglik_unweighted_cpp(occ, bg, as.vector(mu), L)
  wt_val   <- loglik_weighted_local(
    occ, bg, mu, L,
    w1 = rep(1 / n1, n1),
    w2 = rep(1 / n2, n2)
  )

  expect_equal(cpp_val, wt_val, tolerance = 1e-8)
})
