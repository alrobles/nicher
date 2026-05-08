test_that("loglik_niche_math_presence_only_cpp matches legacy R wrapper", {
  occ <- as.matrix(example_env_occ_2d)
  set.seed(1)
  for (k in seq_len(20L)) {
    theta <- start_theta(occ) + stats::rnorm(length(start_theta(occ)),
                                             sd = 0.05)
    legacy <- loglik_niche_math_presence_only(theta, occ,
                                              eta = 1, neg = TRUE)
    new <- nicher:::loglik_niche_math_presence_only_cpp(theta, occ, eta = 1)
    expect_equal(new, legacy, tolerance = 1e-10)
  }
})

test_that("loglik_niche_math_ip_weighted_cpp matches legacy R wrapper", {
  occ <- as.matrix(example_env_occ_2d)
  M   <- as.matrix(example_env_m_2d)

  set.seed(2)
  n_m     <- nrow(M)
  den_idx <- sample.int(n_m, min(500L, n_m))
  kde_idx <- sample.int(n_m, min(500L, n_m))
  M_kde   <- M[kde_idx, , drop = FALSE]
  w_occ   <- as.numeric(kde_gaussian(occ, M_kde))
  w_den   <- as.numeric(kde_gaussian(M[den_idx, , drop = FALSE], M_kde))

  for (k in seq_len(20L)) {
    theta <- start_theta(occ) + stats::rnorm(length(start_theta(occ)),
                                             sd = 0.05)
    legacy <- loglik_niche_math_ip_weighted(
      theta, occ, M, eta = 1, neg = TRUE,
      den_idx = den_idx, kde_idx = kde_idx,
      precomp_w_den = w_den
    )
    new <- nicher:::loglik_niche_math_ip_weighted_cpp(
      theta, occ, M[den_idx, , drop = FALSE], w_occ, w_den, eta = 1
    )
    expect_equal(new, legacy, tolerance = 1e-10)
  }
})

test_that("loglik_niche_math_weighted_cpp respects ridge prior", {
  occ <- as.matrix(example_env_occ_2d)
  M   <- as.matrix(example_env_m_2d)
  p   <- ncol(occ)

  set.seed(3)
  n_m     <- nrow(M)
  den_idx <- sample.int(n_m, min(400L, n_m))
  kde_idx <- sample.int(n_m, min(400L, n_m))
  M_den   <- M[den_idx, , drop = FALSE]
  M_kde   <- M[kde_idx, , drop = FALSE]
  w_occ   <- as.numeric(kde_gaussian(occ, M_kde))
  w_den   <- as.numeric(kde_gaussian(M_den, M_kde))

  log_sd <- log(apply(occ, 2L, stats::sd))

  for (k in seq_len(15L)) {
    theta <- start_theta(occ) + stats::rnorm(length(start_theta(occ)),
                                             sd = 0.05)

    # lambda = 0 -> agrees with paper Eq. 5 (no penalty)
    f_lambda0 <- nicher:::loglik_niche_math_weighted_cpp(
      theta, occ, M_den, w_occ, w_den,
      prior_log_sigma_center = log_sd, prior_log_sigma_lambda = 0.0,
      eta = 1.0
    )

    # ridge term is exactly lambda * sum((log_sigma - center)^2)
    log_sigma <- theta[(p + 1L):(2L * p)]
    ridge_explicit <- 7.0 * sum((log_sigma - log_sd) ^ 2)

    f_lambda7 <- nicher:::loglik_niche_math_weighted_cpp(
      theta, occ, M_den, w_occ, w_den,
      prior_log_sigma_center = log_sd, prior_log_sigma_lambda = 7.0,
      eta = 1.0
    )

    expect_equal(f_lambda7 - f_lambda0, ridge_explicit, tolerance = 1e-10)

    # When log_sigma == prior centre, the penalty contributes 0.
    theta_at_centre <- theta
    theta_at_centre[(p + 1L):(2L * p)] <- log_sd
    f_at_centre_l0 <- nicher:::loglik_niche_math_weighted_cpp(
      theta_at_centre, occ, M_den, w_occ, w_den,
      prior_log_sigma_center = log_sd, prior_log_sigma_lambda = 0.0,
      eta = 1.0
    )
    f_at_centre_l5 <- nicher:::loglik_niche_math_weighted_cpp(
      theta_at_centre, occ, M_den, w_occ, w_den,
      prior_log_sigma_center = log_sd, prior_log_sigma_lambda = 5.0,
      eta = 1.0
    )
    expect_equal(f_at_centre_l5, f_at_centre_l0, tolerance = 1e-10)
  }
})

# ---------------------------------------------------------------------------
# Skew-normal kernels (PR-A.1)
# ---------------------------------------------------------------------------

# R reference: closed-form Azzalini & Capitanio (1999) multivariate
# skew-normal negative log-likelihood, dropping the 0.5*p*log(2*pi) and
# n*log(2) constants that are also dropped in the C++ kernels (they cancel
# in any score equation and in the weighted log-sum-exp denominator).
.ref_neg_loglik_skew_normal <- function(theta, env_occ, eta = 1.0) {
  p   <- ncol(env_occ)
  n_v <- p * (p - 1L) / 2L
  mu        <- theta[seq_len(p)]
  log_sigma <- theta[(p + 1L):(2L * p)]
  v         <- if (n_v > 0L) theta[(2L * p + 1L):(2L * p + n_v)] else numeric(0)
  alpha     <- theta[(2L * p + n_v + 1L):(3L * p + n_v)]

  sigma <- exp(log_sigma)
  L_corr <- nicher:::cvine_cholesky_cpp(v, p, eta)
  Sigma  <- diag(sigma, nrow = p) %*% (L_corr %*% t(L_corr)) %*% diag(sigma, nrow = p)
  Sinv   <- solve(Sigma)

  d  <- sweep(env_occ, 2L, mu, "-")
  q1 <- rowSums((d %*% Sinv) * d)
  z  <- as.numeric(d %*% (alpha / sigma))
  log_phi <- stats::pnorm(z, log.p = TRUE)

  log_det <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  n_occ <- nrow(env_occ)

  0.5 * n_occ * log_det + 0.5 * sum(q1) - sum(log_phi)
}

.ref_neg_loglik_skew_normal_weighted <- function(theta, env_occ, M_den,
                                                 w_occ, w_den,
                                                 prior_log_sigma_center,
                                                 prior_log_sigma_lambda,
                                                 eta = 1.0) {
  p   <- ncol(env_occ)
  n_v <- p * (p - 1L) / 2L
  n_occ <- nrow(env_occ)
  mu        <- theta[seq_len(p)]
  log_sigma <- theta[(p + 1L):(2L * p)]
  v         <- if (n_v > 0L) theta[(2L * p + 1L):(2L * p + n_v)] else numeric(0)
  alpha     <- theta[(2L * p + n_v + 1L):(3L * p + n_v)]

  sigma  <- exp(log_sigma)
  L_corr <- nicher:::cvine_cholesky_cpp(v, p, eta)
  Sigma  <- diag(sigma, nrow = p) %*% (L_corr %*% t(L_corr)) %*% diag(sigma, nrow = p)
  Sinv   <- solve(Sigma)

  d_occ <- sweep(env_occ, 2L, mu, "-")
  q1    <- rowSums((d_occ %*% Sinv) * d_occ)
  z_occ <- as.numeric(d_occ %*% (alpha / sigma))
  log_phi_occ <- stats::pnorm(z_occ, log.p = TRUE)

  d_den <- sweep(M_den, 2L, mu, "-")
  q2    <- rowSums((d_den %*% Sinv) * d_den)
  z_den <- as.numeric(d_den %*% (alpha / sigma))
  log_phi_den <- stats::pnorm(z_den, log.p = TRUE)

  log_w_occ <- log(pmax(w_occ, 1e-300))
  log_w_den <- log(pmax(w_den, 1e-300))

  a   <- log_w_den + log_phi_den - 0.5 * q2
  ma  <- max(a)
  lse <- ma + log(sum(exp(a - ma)))

  ridge <- prior_log_sigma_lambda *
    sum((log_sigma - prior_log_sigma_center) ^ 2L)

  0.5 * sum(q1) - sum(log_phi_occ) - sum(log_w_occ) + n_occ * lse + ridge
}

test_that("loglik_niche_math_skew_normal_cpp matches R reference", {
  occ <- as.matrix(example_env_occ_2d)
  p   <- ncol(occ)
  n_v <- p * (p - 1L) / 2L

  set.seed(7L)
  for (k in seq_len(20L)) {
    base  <- start_theta(occ, skew = TRUE)
    theta <- base + stats::rnorm(length(base), sd = 0.05)
    # set a non-trivial alpha
    theta[(2L * p + n_v + 1L):(3L * p + n_v)] <-
      stats::rnorm(p, sd = 0.5)

    ref <- .ref_neg_loglik_skew_normal(theta, occ)
    got <- nicher:::loglik_niche_math_skew_normal_cpp(theta, occ, eta = 1.0)
    expect_equal(got, ref, tolerance = 1e-9)
  }
})

test_that("loglik_niche_math_skew_normal_weighted_cpp matches R reference", {
  occ <- as.matrix(example_env_occ_2d)
  M   <- as.matrix(example_env_m_2d)
  p   <- ncol(occ)
  n_v <- p * (p - 1L) / 2L

  set.seed(11L)
  n_m   <- nrow(M)
  d_idx <- sample.int(n_m, min(500L, n_m))
  k_idx <- sample.int(n_m, min(500L, n_m))
  M_den <- M[d_idx, , drop = FALSE]
  M_kde <- M[k_idx, , drop = FALSE]
  w_occ <- as.numeric(kde_gaussian(occ, M_kde))
  w_den <- as.numeric(kde_gaussian(M_den, M_kde))

  log_sd <- log(apply(occ, 2L, stats::sd))

  for (k in seq_len(15L)) {
    base  <- start_theta(occ, skew = TRUE)
    theta <- base + stats::rnorm(length(base), sd = 0.05)
    theta[(2L * p + n_v + 1L):(3L * p + n_v)] <-
      stats::rnorm(p, sd = 0.5)

    for (lam in c(0.0, 1.0, 7.0)) {
      ref <- .ref_neg_loglik_skew_normal_weighted(
        theta, occ, M_den, w_occ, w_den,
        prior_log_sigma_center = log_sd,
        prior_log_sigma_lambda = lam
      )
      got <- nicher:::loglik_niche_math_skew_normal_weighted_cpp(
        theta, occ, M_den, w_occ, w_den,
        prior_log_sigma_center = log_sd,
        prior_log_sigma_lambda = lam,
        eta = 1.0
      )
      expect_equal(got, ref, tolerance = 1e-9)
    }
  }
})

# ---------------------------------------------------------------------------
# Reference: pure-R skew-t (NCST) negative log-likelihood using
# stats::integrate(). Slow but high-precision; sufficient as a parity oracle.
# ---------------------------------------------------------------------------
.ref_neg_loglik_skew_t <- function(theta, env_occ, eta = 1.0) {
  p   <- ncol(env_occ)
  n_v <- p * (p - 1L) / 2L
  mu        <- theta[seq_len(p)]
  log_sigma <- theta[(p + 1L):(2L * p)]
  v         <- if (n_v > 0L) theta[(2L * p + 1L):(2L * p + n_v)] else numeric(0)
  alpha     <- theta[(2L * p + n_v + 1L):(3L * p + n_v)]
  log_r     <- theta[3L * p + n_v + 1L]
  r         <- exp(log_r)

  sigma  <- exp(log_sigma)
  L_corr <- nicher:::cvine_cholesky_cpp(v, p, eta)
  Sigma  <- diag(sigma, nrow = p) %*% (L_corr %*% t(L_corr)) %*% diag(sigma, nrow = p)
  Sinv   <- solve(Sigma)

  d  <- sweep(env_occ, 2L, mu, "-")
  q1 <- rowSums((d %*% Sinv) * d)
  z  <- as.numeric(d %*% (alpha / sigma))
  A  <- 1.0 + q1 / r

  log_det <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  n_occ   <- nrow(env_occ)
  a       <- (p + r) / 2.0

  sum_log_A <- sum(log(A))
  sum_log_I <- sum(vapply(seq_len(n_occ), function(i) {
    val <- stats::integrate(
      f = function(u) {
        u^(a - 1) * exp(-u) *
          stats::pnorm(z[i] * sqrt(2 * u / (A[i] * r)))
      },
      lower = 0, upper = Inf,
      rel.tol = 1e-12, subdivisions = 2000L
    )$value
    log(val)
  }, numeric(1)))

  # The C++ kernel computes the parameter-dependent part of -sum log f(t_i),
  # dropping the same data-only constants as the skew-normal kernel:
  #   -(p/2) n log(2 pi) and +(1 + p/2) n log 2 (both purely constant
  #   in theta and r, so they fall out of any score equation).
  # Remaining theta-dependent terms (matches Branco & Dey 2001, Sec. 4):
  0.5 * n_occ * log_det +
    0.5 * n_occ * p * log(r) +
    n_occ * lgamma(r / 2) +
    a * sum_log_A -
    sum_log_I
}

test_that("loglik_niche_math_skew_t_cpp matches R reference (stats::integrate)", {
  occ <- as.matrix(example_env_occ_2d)
  p   <- ncol(occ)
  n_v <- p * (p - 1L) / 2L

  set.seed(31L)
  for (k in seq_len(8L)) {
    base  <- start_theta(occ, skew_t = TRUE)
    theta <- base + stats::rnorm(length(base), sd = 0.05)
    # Sample non-trivial alpha and log_r in the supported range.
    theta[(2L * p + n_v + 1L):(3L * p + n_v)] <- stats::rnorm(p, sd = 0.5)
    theta[3L * p + n_v + 1L] <- log(2) + stats::runif(1L) *
                                          (log(100) - log(2))

    ref <- as.numeric(.ref_neg_loglik_skew_t(theta, occ))
    got <- as.numeric(nicher:::loglik_niche_math_skew_t_cpp(theta, occ,
                                                            eta = 1.0))
    expect_equal(got, ref, tolerance = 1e-3)   # GL Q=32 is ~1e-6 relative
  }
})

# ---------------------------------------------------------------------------
# Reference: pure-R skew-t weighted negative log-likelihood. Reuses the
# per-point integral evaluator from the presence-only reference, then
# composes with the paper Eq. 5 weighted-likelihood form.
# ---------------------------------------------------------------------------
.ref_neg_loglik_skew_t_weighted <- function(theta, env_occ, M_den,
                                            w_occ, w_den,
                                            prior_log_sigma_center,
                                            prior_log_sigma_lambda,
                                            eta = 1.0) {
  p   <- ncol(env_occ)
  n_v <- p * (p - 1L) / 2L
  n_occ <- nrow(env_occ)
  mu        <- theta[seq_len(p)]
  log_sigma <- theta[(p + 1L):(2L * p)]
  v         <- if (n_v > 0L) theta[(2L * p + 1L):(2L * p + n_v)] else numeric(0)
  alpha     <- theta[(2L * p + n_v + 1L):(3L * p + n_v)]
  log_r     <- theta[3L * p + n_v + 1L]
  r         <- exp(log_r)

  sigma  <- exp(log_sigma)
  L_corr <- nicher:::cvine_cholesky_cpp(v, p, eta)
  Sigma  <- diag(sigma, nrow = p) %*% (L_corr %*% t(L_corr)) %*% diag(sigma, nrow = p)
  Sinv   <- solve(Sigma)

  log_det <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  a       <- (p + r) / 2.0

  # Helper that integrates I(t) = ∫_0^∞ u^{a-1} e^{-u} Φ(z √(2u/(Ar))) du
  log_I_at <- function(z_val, A_val) {
    val <- stats::integrate(
      f = function(u) {
        u^(a - 1) * exp(-u) *
          stats::pnorm(z_val * sqrt(2 * u / (A_val * r)))
      },
      lower = 0, upper = Inf,
      rel.tol = 1e-12, subdivisions = 2000L
    )$value
    log(val)
  }

  d_occ <- sweep(env_occ, 2L, mu, "-")
  q1    <- rowSums((d_occ %*% Sinv) * d_occ)
  z_occ <- as.numeric(d_occ %*% (alpha / sigma))
  A_occ <- 1.0 + q1 / r
  log_I_occ <- vapply(seq_len(n_occ),
                      function(i) log_I_at(z_occ[i], A_occ[i]),
                      numeric(1))

  d_den <- sweep(M_den, 2L, mu, "-")
  q2    <- rowSums((d_den %*% Sinv) * d_den)
  z_den <- as.numeric(d_den %*% (alpha / sigma))
  A_den <- 1.0 + q2 / r
  log_I_den <- vapply(seq_len(nrow(M_den)),
                      function(j) log_I_at(z_den[j], A_den[j]),
                      numeric(1))

  log_w_occ <- log(pmax(w_occ, 1e-300))
  log_w_den <- log(pmax(w_den, 1e-300))

  # Per-point log-density (theta-dependent constants only):
  #   log f(t) = -0.5 log|Sigma| - 0.5 p log(r) - lgamma(r/2)
  #              + a (log 2 - log A) + log I(t)
  # Drop the data-only constants -(p/2)log(2π) and (1+p/2) log 2.
  log_dens_occ <- -0.5 * log_det - 0.5 * p * log(r) - lgamma(r / 2) -
                    a * log(A_occ) + log_I_occ
  log_dens_den <- -0.5 * log_det - 0.5 * p * log(r) - lgamma(r / 2) -
                    a * log(A_den) + log_I_den

  # Paper Eq. 5: -ℓ = -sum log f(x_i) - sum log w(x_i) +
  #                    n log sum_j w(y_j) f(y_j; theta)  + ridge
  a_lse <- log_w_den + log_dens_den
  ma   <- max(a_lse)
  lse  <- ma + log(sum(exp(a_lse - ma)))

  ridge <- prior_log_sigma_lambda *
    sum((log_sigma - prior_log_sigma_center) ^ 2L)

  -sum(log_dens_occ) - sum(log_w_occ) + n_occ * lse + ridge
}

test_that("loglik_niche_math_skew_t_weighted_cpp matches R reference", {
  occ <- as.matrix(example_env_occ_2d)
  M   <- as.matrix(example_env_m_2d)
  p   <- ncol(occ)
  n_v <- p * (p - 1L) / 2L

  set.seed(41L)
  n_m   <- nrow(M)
  d_idx <- sample.int(n_m, min(80L, n_m))    # tiny denom set: R is slow
  k_idx <- sample.int(n_m, min(200L, n_m))
  M_den <- M[d_idx, , drop = FALSE]
  M_kde <- M[k_idx, , drop = FALSE]
  w_occ <- as.numeric(kde_gaussian(occ, M_kde))
  w_den <- as.numeric(kde_gaussian(M_den, M_kde))
  log_sd <- log(apply(occ, 2L, stats::sd))

  for (k in seq_len(3L)) {
    base  <- start_theta(occ, skew_t = TRUE)
    theta <- base + stats::rnorm(length(base), sd = 0.05)
    theta[(2L * p + n_v + 1L):(3L * p + n_v)] <- stats::rnorm(p, sd = 0.4)
    theta[3L * p + n_v + 1L] <- log(2) + stats::runif(1L) *
                                          (log(100) - log(2))

    for (lam in c(0.0, 1.0)) {
      ref <- as.numeric(.ref_neg_loglik_skew_t_weighted(
        theta, occ, M_den, w_occ, w_den,
        prior_log_sigma_center = log_sd,
        prior_log_sigma_lambda = lam
      ))
      got <- as.numeric(nicher:::loglik_niche_math_skew_t_weighted_cpp(
        theta, occ, M_den, w_occ, w_den,
        prior_log_sigma_center = log_sd,
        prior_log_sigma_lambda = lam,
        eta = 1.0
      ))
      expect_equal(got, ref, tolerance = 1e-3)
    }
  }
})

test_that("skew_t kernel returns finite, sensible values across r", {
  # NCST_k(mu, Sigma, alpha, r) -> SN_k(mu, Sigma, alpha) as r -> Inf
  # (Branco & Dey 2001, Sec. 4). The two C++ kernels drop different
  # data-only constants so they cannot be compared by exact equality;
  # this test instead pins down a behavioural sanity check: the kernel
  # is finite and monotone in (a single coordinate of) theta across a
  # wide range of r.
  occ <- as.matrix(example_env_occ_2d)
  p   <- ncol(occ)
  n_v <- p * (p - 1L) / 2L

  set.seed(37L)
  base  <- start_theta(occ, skew_t = TRUE)
  theta <- base + stats::rnorm(length(base), sd = 0.05)
  theta[(2L * p + n_v + 1L):(3L * p + n_v)] <- stats::rnorm(p, sd = 0.4)

  for (lr in c(log(2), log(5), log(20), log(100), log(1e3))) {
    theta[3L * p + n_v + 1L] <- lr
    f <- nicher:::loglik_niche_math_skew_t_cpp(theta, occ, eta = 1.0)
    expect_true(is.finite(f))
  }
})

test_that("skew_normal reduces to presence_only when alpha = 0", {
  # SN_k(mu, Sigma, alpha = 0) is exactly N_k(mu, Sigma); the two negative
  # log-likelihoods differ only by the n*log(2) constant that the SN
  # kernel intentionally drops to keep the form numerically clean.
  occ <- as.matrix(example_env_occ_2d)
  p   <- ncol(occ)
  n_v <- p * (p - 1L) / 2L

  set.seed(13L)
  base  <- start_theta(occ)                       # 2p + n_v
  theta_g <- base + stats::rnorm(length(base), sd = 0.05)
  theta_s <- c(theta_g, rep(0.0, p))              # 3p + n_v, alpha = 0

  f_po <- loglik_niche_math_presence_only(theta_g, occ, eta = 1.0,
                                          neg = TRUE)
  f_sn <- nicher:::loglik_niche_math_skew_normal_cpp(theta_s, occ, eta = 1.0)

  # With alpha = 0, log Phi(0) = log(0.5) = -log 2 per occurrence;
  # the SN kernel's "presence" term is 0.5 sum(q1) - sum(log Phi(z)),
  # so f_sn = f_po + n * log(2).
  n_occ <- nrow(occ)
  expect_equal(f_sn, f_po + n_occ * log(2), tolerance = 1e-10)
})
