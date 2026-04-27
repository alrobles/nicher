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

test_that("loglik_niche_math_weighted_cpp matches legacy R wrapper", {
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
    legacy <- loglik_niche_math_weighted(
      theta, occ, M, eta = 1, neg = TRUE,
      den_idx = den_idx, kde_idx = kde_idx,
      precomp_w_den = w_den
    )
    new <- nicher:::loglik_niche_math_weighted_cpp(
      theta, occ, M[den_idx, , drop = FALSE], w_occ, w_den, eta = 1
    )
    expect_equal(new, legacy, tolerance = 1e-10)
  }
})
