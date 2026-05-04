# Cross-validation: cv_nicher() formula sanity + fingerprint validation.

.fit_w_small <- function(seed = 1L) {
  set.seed(seed)
  optimize_niche(
    env_occ = example_env_occ_2d, env_m = example_env_m_2d,
    num_starts = 5L, breadth = 0.1,
    likelihood = "weighted",
    prior_mu_lambda = 1, prior_log_sigma_lambda = 1, prior_alpha_lambda = 0,
    verbose = FALSE
  )
}
.fit_po_small <- function(seed = 1L) {
  set.seed(seed)
  optimize_niche(
    env_occ = example_env_occ_2d, env_m = NULL,
    num_starts = 5L, breadth = 0.1,
    likelihood = "presence_only", verbose = FALSE
  )
}

# -----------------------------------------------------------------------
test_that("cv_nicher returns a sane nicher_cv list (kfold, weighted)", {
  fit <- .fit_w_small()
  cv  <- cv_nicher(fit, example_env_occ_2d, example_env_m_2d,
                   type = "kfold", k = 5L, seed = 42L)
  expect_s3_class(cv, "nicher_cv")
  expect_true(is.finite(cv$cv_loglik))
  expect_equal(cv$type, "kfold")
  expect_equal(cv$k,    5L)
  expect_equal(nrow(cv$per_fold), 5L)
  expect_equal(sum(cv$per_fold$n_test), nrow(example_env_occ_2d))
  expect_equal(cv$cv_loglik_mean, cv$cv_loglik / nrow(example_env_occ_2d),
               tolerance = 1e-12)
})

test_that("cv_nicher works for presence_only (no env_m)", {
  fit <- .fit_po_small()
  cv  <- cv_nicher(fit, example_env_occ_2d, env_m = NULL,
                   type = "kfold", k = 5L, seed = 42L)
  expect_s3_class(cv, "nicher_cv")
  expect_true(is.finite(cv$cv_loglik))
  expect_equal(cv$type, "kfold")
})

test_that("cv_nicher refuses mismatched env_occ", {
  fit <- .fit_w_small()
  bad <- example_env_occ_2d
  bad[1L, 1L] <- bad[1L, 1L] + 100
  expect_error(
    cv_nicher(fit, bad, example_env_m_2d, type = "kfold", k = 5L),
    regexp = "fingerprint mismatch"
  )
})

test_that("cv_nicher refuses mismatched env_m", {
  fit <- .fit_w_small()
  bad_m <- example_env_m_2d
  bad_m[1L, 1L] <- bad_m[1L, 1L] + 100
  expect_error(
    cv_nicher(fit, example_env_occ_2d, bad_m, type = "kfold", k = 5L),
    regexp = "fingerprint mismatch"
  )
})

test_that("cv_nicher refuses missing env_m for weighted family", {
  fit <- .fit_w_small()
  expect_error(
    cv_nicher(fit, example_env_occ_2d, env_m = NULL, type = "kfold", k = 5L),
    regexp = "env_m.*required"
  )
})

test_that("cv_nicher LOO basic structure (small data)", {
  fit <- .fit_po_small()
  small_occ <- example_env_occ_2d[seq_len(20L), , drop = FALSE]
  fit_small <- optimize_niche(
    env_occ = small_occ, env_m = NULL,
    num_starts = 3L, breadth = 0.1,
    likelihood = "presence_only", verbose = FALSE
  )
  cv <- cv_nicher(fit_small, small_occ, env_m = NULL, type = "loo")
  expect_s3_class(cv, "nicher_cv")
  expect_equal(cv$type, "loo")
  expect_equal(cv$k, 20L)
  expect_equal(nrow(cv$per_fold), 20L)
  expect_true(all(cv$per_fold$n_test == 1L))
})

test_that("cv_nicher errors clearly on pre-3.3 nicher object", {
  fit <- .fit_po_small()
  fit$fit_args <- NULL                              # simulate old object
  expect_error(
    cv_nicher(fit, example_env_occ_2d, env_m = NULL, type = "kfold", k = 5L),
    regexp = "lacks `fit_args`"
  )
})

test_that("compare_nicher accepts comparison_basis = 'penalised'", {
  fit_a <- .fit_w_small(seed = 1L)
  set.seed(2L)
  fit_b <- optimize_niche(
    env_occ = example_env_occ_2d, env_m = example_env_m_2d,
    num_starts = 5L, breadth = 0.1, likelihood = "weighted",
    prior_mu_lambda = 5, prior_log_sigma_lambda = 5,
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
