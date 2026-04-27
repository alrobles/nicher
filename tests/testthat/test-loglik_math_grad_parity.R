test_that("hybrid analytic gradient matches central FD", {
  occ <- as.matrix(example_env_occ_2d)
  M   <- as.matrix(example_env_m_2d)

  set.seed(3)
  n_m     <- nrow(M)
  den_idx <- sample.int(n_m, min(500L, n_m))
  kde_idx <- sample.int(n_m, min(500L, n_m))
  M_kde   <- M[kde_idx, , drop = FALSE]
  w_occ   <- as.numeric(kde_gaussian(occ, M_kde))
  w_den   <- as.numeric(kde_gaussian(M[den_idx, , drop = FALSE], M_kde))
  M_den   <- M[den_idx, , drop = FALSE]

  fd_grad <- function(theta) {
    n <- length(theta)
    g <- numeric(n)
    for (i in seq_len(n)) {
      dx <- abs(theta[i]) * 1e-6 + 1e-8
      th_p <- theta; th_p[i] <- theta[i] + dx
      th_m <- theta; th_m[i] <- theta[i] - dx
      fp <- nicher:::loglik_niche_math_weighted_cpp(
        th_p, occ, M_den, w_occ, w_den, eta = 1)
      fm <- nicher:::loglik_niche_math_weighted_cpp(
        th_m, occ, M_den, w_occ, w_den, eta = 1)
      g[i] <- (fp - fm) / (2 * dx)
    }
    g
  }

  for (k in seq_len(5L)) {
    theta <- start_theta(occ) + stats::rnorm(length(start_theta(occ)),
                                             sd = 0.05)
    res <- nicher:::loglik_niche_math_weighted_grad_cpp(
      theta, occ, M_den, w_occ, w_den, eta = 1)
    g_fd <- fd_grad(theta)

    # Tolerance dominated by FD precision (~1e-6 relative)
    expect_equal(res$gradient, g_fd, tolerance = 1e-5)
  }
})
