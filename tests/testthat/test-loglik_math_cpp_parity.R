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

test_that("loglik_niche_math_kde_bias_corrected_cpp matches legacy R wrapper", {
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
    legacy <- loglik_niche_math_kde_bias_corrected(
      theta, occ, M, eta = 1, neg = TRUE,
      den_idx = den_idx, kde_idx = kde_idx,
      precomp_w_den = w_den
    )
    new <- nicher:::loglik_niche_math_kde_bias_corrected_cpp(
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
