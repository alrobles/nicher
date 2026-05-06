test_that("penalised-MLE defaults keep weighted mu near the PO mu", {
  skip_if_not_installed("pomp")

  env_occ <- example_env_occ_2d
  env_m   <- example_env_m_2d
  stopifnot(nrow(env_occ) > 0, nrow(env_m) > 0)

  po <- optimize_niche(
    env_occ    = env_occ,
    env_m      = NULL,
    num_starts = 3L,
    breadth    = 0.1,
    likelihood = "presence_only",
    seed       = 42L,
    control    = list(maxeval = 120L),
    verbose    = FALSE
  )
  expect_true(po$best$convergence %in% c(1L, 2L))
  p <- ncol(env_occ)
  mu_po    <- as.numeric(po$best$theta[seq_len(p)])
  sigma_po <- exp(as.numeric(po$best$theta[p + seq_len(p)]))

  w <- optimize_niche(
    env_occ    = env_occ,
    env_m      = env_m,
    num_starts = 3L,
    breadth    = 0.1,
    likelihood = "weighted",
    seed       = 42L,
    control    = list(maxeval = 120L),
    verbose    = FALSE
  )
  expect_true(w$best$convergence %in% c(1L, 2L))
  mu_w <- as.numeric(w$best$theta[seq_len(p)])

  # Drift in PO-sigma units. The regression bound is 5 to give some
  # optimizer slack.
  drift <- abs(mu_w - mu_po) / sigma_po
  expect_true(all(drift < 5),
              info = sprintf("PO-sigma drift: %s", paste(round(drift, 3),
                                                        collapse = ", ")))
})

test_that("prior_alpha_lambda bounds alpha in presence-only skew_normal", {
  skip_if_not_installed("pomp")
  env_occ <- example_env_occ_2d

  fit <- optimize_niche(
    env_occ    = env_occ,
    env_m      = NULL,
    num_starts = 3L,
    breadth    = 0.1,
    likelihood = "skew_normal",
    seed       = 42L,
    control    = list(maxeval = 120L),
    verbose    = FALSE
  )
  expect_true(fit$best$convergence %in% c(1L, 2L))
  p <- ncol(env_occ)
  # theta = [mu(p), log_sigma(p), v(p(p-1)/2), alpha(p)]
  n_v <- p * (p - 1L) / 2L
  alpha <- as.numeric(fit$best$theta[2L * p + n_v + seq_len(p)])
  expect_true(all(abs(alpha) < 100),
              info = sprintf("alpha: %s", paste(round(alpha, 3),
                                                collapse = ", ")))
})
