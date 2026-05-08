# Cross-validation: cv_nicher() formula sanity + fingerprint validation.

.fit_w_small <- function(seed = 1L) {
  set.seed(seed)
  optimize_niche(
    env_occ = example_env_occ_2d[seq_len(24L), , drop = FALSE],
    env_m = example_env_m_2d,
    num_starts = 2L, breadth = 0.1,
    likelihood = "weighted",
    prior_mu_lambda = 1, prior_log_sigma_lambda = 1, prior_alpha_lambda = 0,
    m_subsample = 500L, m_kde_subsample = 500L,
    control = list(maxeval = 80L),
    verbose = FALSE
  )
}
.fit_po_small <- function(seed = 1L) {
  set.seed(seed)
  optimize_niche(
    env_occ = example_env_occ_2d[seq_len(24L), , drop = FALSE],
    env_m = NULL,
    num_starts = 2L, breadth = 0.1,
    likelihood = "presence_only",
    control = list(maxeval = 80L),
    verbose = FALSE
  )
}

# -----------------------------------------------------------------------
test_that("cv_nicher returns a sane nicher_cv list (kfold, weighted)", {
  fit <- .fit_w_small()
  occ <- example_env_occ_2d[seq_len(24L), , drop = FALSE]
  M <- example_env_m_2d
  cv  <- cv_nicher(fit, occ, M, type = "kfold", k = 3L, seed = 42L,
                   ucminf_control = list(maxeval = 80L))
  expect_s3_class(cv, "nicher_cv")
  expect_true(is.finite(cv$cv_loglik))
  expect_equal(cv$type, "kfold")
  expect_equal(cv$k,    3L)
  expect_equal(nrow(cv$per_fold), 3L)
  expect_equal(sum(cv$per_fold$n_test), nrow(occ))
  expect_equal(cv$cv_loglik_mean, cv$cv_loglik / nrow(occ),
               tolerance = 1e-12)
})

test_that("cv_nicher works for presence_only (no env_m)", {
  fit <- .fit_po_small()
  occ <- example_env_occ_2d[seq_len(24L), , drop = FALSE]
  cv  <- cv_nicher(fit, occ, env_m = NULL, type = "kfold", k = 3L,
                   seed = 42L, ucminf_control = list(maxeval = 80L))
  expect_s3_class(cv, "nicher_cv")
  expect_true(is.finite(cv$cv_loglik))
  expect_equal(cv$type, "kfold")
})

test_that("cv_nicher refuses mismatched env_occ", {
  fit <- .fit_w_small()
  bad <- example_env_occ_2d[seq_len(24L), , drop = FALSE]
  M <- example_env_m_2d
  bad[1L, 1L] <- bad[1L, 1L] + 100
  expect_error(
    cv_nicher(fit, bad, M, type = "kfold", k = 3L),
    regexp = "fingerprint mismatch"
  )
})

test_that("cv_nicher refuses mismatched env_m", {
  fit <- .fit_w_small()
  occ <- example_env_occ_2d[seq_len(24L), , drop = FALSE]
  bad_m <- example_env_m_2d
  bad_m[1L, 1L] <- bad_m[1L, 1L] + 100
  expect_error(
    cv_nicher(fit, occ, bad_m, type = "kfold", k = 3L),
    regexp = "fingerprint mismatch"
  )
})

test_that("cv_nicher refuses missing env_m for weighted family", {
  fit <- .fit_w_small()
  expect_error(
    cv_nicher(fit, example_env_occ_2d[seq_len(24L), , drop = FALSE],
              env_m = NULL, type = "kfold", k = 3L),
    regexp = "env_m.*required"
  )
})

test_that("cv_nicher LOO basic structure (tiny data)", {
  fit <- .fit_po_small()
  small_occ <- example_env_occ_2d[seq_len(6L), , drop = FALSE]
  fit_small <- optimize_niche(
    env_occ = small_occ, env_m = NULL,
    num_starts = 2L, breadth = 0.1,
    likelihood = "presence_only", control = list(maxeval = 60L),
    verbose = FALSE
  )
  cv <- cv_nicher(fit_small, small_occ, env_m = NULL, type = "loo",
                  ucminf_control = list(maxeval = 40L))
  expect_s3_class(cv, "nicher_cv")
  expect_equal(cv$type, "loo")
  expect_equal(cv$k, 6L)
  expect_equal(nrow(cv$per_fold), 6L)
  expect_true(all(cv$per_fold$n_test == 1L))
})

test_that("cv_nicher errors clearly on pre-3.3 nicher object", {
  fit <- .fit_po_small()
  fit$fit_args <- NULL                              # simulate old object
  expect_error(
    cv_nicher(fit, example_env_occ_2d[seq_len(24L), , drop = FALSE],
              env_m = NULL, type = "kfold", k = 3L),
    regexp = "lacks `fit_args`"
  )
})

test_that("compare_nicher accepts comparison_basis = 'penalised'", {
  fit_a <- .fit_w_small(seed = 1L)
  set.seed(2L)
  fit_b <- optimize_niche(
    env_occ = example_env_occ_2d[seq_len(24L), , drop = FALSE],
    env_m = example_env_m_2d,
    num_starts = 2L, breadth = 0.1, likelihood = "weighted",
    prior_mu_lambda = 5, prior_log_sigma_lambda = 5,
    m_subsample = 500L, m_kde_subsample = 500L,
    control = list(maxeval = 80L),
    verbose = FALSE
  )
  cmp_un  <- compare_nicher(a = fit_a, b = fit_b, comparison_basis = "unpenalised")
  cmp_pen <- compare_nicher(a = fit_a, b = fit_b, comparison_basis = "penalised")
  expect_s3_class(cmp_un,  "data.frame")
  expect_s3_class(cmp_pen, "data.frame")
  expect_identical(attr(cmp_un,  "comparison_basis"), "unpenalised")
  expect_identical(attr(cmp_pen, "comparison_basis"), "penalised")
  # Penalised AIC differs from unpenalised AIC for penalised fits.
  # (It would be identical only if both fits had penalty == 0.)
  expect_false(isTRUE(all.equal(cmp_un$AIC, cmp_pen$AIC)))
})
