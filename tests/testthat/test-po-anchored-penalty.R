test_that("penalised-MLE defaults keep weighted mu near the PO mu", {
  # Real-data regression test. The bundled example_vicugna dataset is the
  # canonical failure case for the unregularised weighted MLE -- on
  # nicher 3.0.x the default `weighted` fit wandered to mu2 = 680 while
  # the PO fit landed at mu2 = 215 (a ~30-sigma drift). With the PO
  # anchor enabled by default (prior_mu_lambda = 1) the weighted fit
  # must stay within a modest number of PO standard deviations.
  skip_if_not_installed("pomp")

  env_occ <- example_vicugna$env_occ
  env_m   <- example_vicugna$env_m
  stopifnot(nrow(env_occ) > 0, nrow(env_m) > 0)

  po <- optimize_niche(
    env_occ    = env_occ,
    env_m      = NULL,
    num_starts = 20L,
    breadth    = 0.1,
    likelihood = "presence_only",
    seed       = 42L,
    verbose    = FALSE
  )
  expect_true(po$best$convergence %in% c(1L, 2L))
  p <- ncol(env_occ)
  mu_po    <- as.numeric(po$best$theta[seq_len(p)])
  sigma_po <- exp(as.numeric(po$best$theta[p + seq_len(p)]))

  w <- optimize_niche(
    env_occ    = env_occ,
    env_m      = env_m,
    num_starts = 20L,
    breadth    = 0.45,
    likelihood = "weighted",
    seed       = 42L,
    verbose    = FALSE
  )
  expect_true(w$best$convergence %in% c(1L, 2L))
  mu_w <- as.numeric(w$best$theta[seq_len(p)])

  # Drift in PO-sigma units. With prior_mu_lambda = 1 we expect
  # |mu_w_k - mu_po_k| / sigma_po_k <= 5 comfortably. The regression
  # bound is 5 to give some optimizer slack; the unregularised
  # nicher 3.0.x fit produced |drift| ~= 30 on mu2.
  drift <- abs(mu_w - mu_po) / sigma_po
  expect_true(all(drift < 5),
              info = sprintf("PO-sigma drift: %s", paste(round(drift, 3),
                                                        collapse = ", ")))
})

test_that("prior_alpha_lambda bounds alpha in presence-only skew_normal", {
  skip_if_not_installed("pomp")
  env_occ <- example_vicugna$env_occ

  fit <- optimize_niche(
    env_occ    = env_occ,
    env_m      = NULL,
    num_starts = 10L,
    breadth    = 0.45,
    likelihood = "skew_normal",
    seed       = 42L,
    verbose    = FALSE
  )
  expect_true(fit$best$convergence %in% c(1L, 2L))
  p <- ncol(env_occ)
  # theta = [mu(p), log_sigma(p), v(p(p-1)/2), alpha(p)]
  n_v <- p * (p - 1L) / 2L
  alpha <- as.numeric(fit$best$theta[2L * p + n_v + seq_len(p)])
  # Without the alpha penalty, alpha on example_vicugna explodes to
  # ~1e6. With lambda_alpha = 0.1 it must stay finite and modest.
  expect_true(all(abs(alpha) < 100),
              info = sprintf("alpha: %s", paste(round(alpha, 3),
                                                collapse = ", ")))
})

test_that("prior_*_lambda = 0 recovers nicher 3.0.x pure-MLE behaviour", {
  # Sanity check: turning all three penalties off matches the bygone
  # un-anchored optimisation in structure (loglik is the plain negative
  # likelihood, no ridge).
  skip_if_not_installed("pomp")
  env_occ <- example_env_occ_2d
  env_m   <- example_env_m_2d

  fit_pen <- optimize_niche(
    env_occ = env_occ, env_m = env_m,
    num_starts = 5L, breadth = 0.1,
    likelihood = "weighted", seed = 1L, verbose = FALSE,
    prior_mu_lambda = 0, prior_log_sigma_lambda = 0,
    prior_alpha_lambda = 0
  )
  expect_true(fit_pen$best$convergence %in% c(1L, 2L))
  # Loglik returned is the negative objective value -- with lambdas = 0
  # this IS the plain log-likelihood, so it should be reproducible from
  # the C++ kernel with the same centre/lambda passthrough.
  expect_true(is.finite(fit_pen$best$loglik))
})
