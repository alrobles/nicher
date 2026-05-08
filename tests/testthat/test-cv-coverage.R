# Additional coverage tests for cv.R and ic.R

# Shared lightweight fit helper
.cv_cov_fit_po <- function(n = 20L, seed = 10L) {
  set.seed(seed)
  optimize_niche(
    env_occ = example_env_occ_2d[seq_len(n), , drop = FALSE],
    env_m = NULL, num_starts = 2L, breadth = 0.1,
    likelihood = "presence_only",
    control = list(maxeval = 60L), verbose = FALSE
  )
}

.cv_cov_fit_w <- function(n = 20L, seed = 10L) {
  set.seed(seed)
  optimize_niche(
    env_occ = example_env_occ_2d[seq_len(n), , drop = FALSE],
    env_m = example_env_m_2d, num_starts = 2L, breadth = 0.1,
    likelihood = "weighted",
    prior_mu_lambda = 1, prior_log_sigma_lambda = 1,
    m_subsample = 500L, m_kde_subsample = 500L,
    control = list(maxeval = 60L), verbose = FALSE
  )
}

# ---- cv_nicher input validation ------------------------------------------

test_that("cv_nicher rejects non-nicher input", {
  expect_error(cv_nicher("not_nicher", example_env_occ_2d),
               "must be a `nicher` object")
})

test_that("cv_nicher rejects bad k values", {
  fit <- .cv_cov_fit_po()
  occ <- example_env_occ_2d[seq_len(20L), , drop = FALSE]
  expect_error(cv_nicher(fit, occ, type = "kfold", k = 1L), "must be an integer")
  expect_error(cv_nicher(fit, occ, type = "kfold", k = 100L), "must be an integer")
})

test_that("cv_nicher rejects bad num_starts_cv", {
  fit <- .cv_cov_fit_po()
  occ <- example_env_occ_2d[seq_len(20L), , drop = FALSE]
  expect_error(cv_nicher(fit, occ, type = "kfold", k = 3L, num_starts_cv = -1),
               "num_starts_cv.*positive")
})

test_that("cv_nicher multi-start fold refit works", {
  fit <- .cv_cov_fit_po(n = 12L, seed = 20L)
  occ <- example_env_occ_2d[seq_len(12L), , drop = FALSE]
  cv <- cv_nicher(fit, occ, env_m = NULL, type = "kfold", k = 3L,
                  seed = 42L, num_starts_cv = 2L,
                  ucminf_control = list(maxeval = 40L))
  expect_s3_class(cv, "nicher_cv")
  expect_true(is.finite(cv$cv_loglik))
})

# ---- print.nicher_cv -----------------------------------------------------

test_that("print.nicher_cv returns invisibly without error", {
  fit <- .cv_cov_fit_po()
  occ <- example_env_occ_2d[seq_len(20L), , drop = FALSE]
  cv <- cv_nicher(fit, occ, env_m = NULL, type = "kfold", k = 3L,
                  seed = 42L, ucminf_control = list(maxeval = 60L))
  out <- capture.output(ret <- print(cv))
  expect_identical(ret, cv)
  expect_true(any(grepl("nicher_cv", out)))
})

# ---- logLik.nicher / nobs.nicher -----------------------------------------

test_that("logLik.nicher returns correct structure", {
  fit <- .cv_cov_fit_po()
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_true(is.finite(as.numeric(ll)))
  expect_true(!is.null(attr(ll, "df")))
  expect_true(!is.null(attr(ll, "nobs")))
  expect_equal(attr(ll, "nobs"), 20L)
})

test_that("logLik.nicher errors on old object", {
  fit <- .cv_cov_fit_po()
  fit$best$loglik_unpenalised <- NULL
  expect_error(logLik(fit), "loglik_unpenalised")
})

test_that("nobs.nicher returns correct count", {
  fit <- .cv_cov_fit_po()
  expect_equal(nobs(fit), 20L)
})

test_that("nobs.nicher errors on old object", {
  fit <- .cv_cov_fit_po()
  fit$best$nobs <- NULL
  expect_error(nobs(fit), "nobs")
})

# ---- AIC / BIC -----------------------------------------------------------

test_that("AIC.nicher and BIC.nicher return numeric for single fit", {
  fit <- .cv_cov_fit_po()
  expect_true(is.numeric(AIC(fit)))
  expect_true(is.numeric(BIC(fit)))
  expect_true(is.finite(AIC(fit)))
  expect_true(is.finite(BIC(fit)))
})

test_that("AIC.nicher returns a data.frame for multiple fits", {
  fit_a <- .cv_cov_fit_po(seed = 1L)
  fit_b <- .cv_cov_fit_po(seed = 2L)
  out <- AIC(fit_a, fit_b)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2L)
  expect_true("AIC" %in% names(out))
})

test_that("BIC.nicher returns a data.frame for multiple fits", {
  fit_a <- .cv_cov_fit_po(seed = 1L)
  fit_b <- .cv_cov_fit_po(seed = 2L)
  out <- BIC(fit_a, fit_b)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2L)
  expect_true("BIC" %in% names(out))
})

# ---- compare_nicher ------------------------------------------------------

test_that("compare_nicher works with named list input", {
  fit_a <- .cv_cov_fit_po(seed = 1L)
  fit_b <- .cv_cov_fit_po(seed = 2L)
  cmp <- compare_nicher(list(a = fit_a, b = fit_b))
  expect_s3_class(cmp, "data.frame")
  expect_equal(nrow(cmp), 2L)
})

test_that("compare_nicher supports sort_by = 'BIC' and 'loglik'", {
  fit_a <- .cv_cov_fit_po(seed = 1L)
  fit_b <- .cv_cov_fit_po(seed = 2L)
  cmp_bic <- compare_nicher(fit_a, fit_b, sort_by = "BIC")
  cmp_ll  <- compare_nicher(fit_a, fit_b, sort_by = "loglik")
  expect_s3_class(cmp_bic, "data.frame")
  expect_s3_class(cmp_ll, "data.frame")
})

test_that("compare_nicher rejects < 2 fits", {
  fit <- .cv_cov_fit_po()
  expect_error(compare_nicher(fit), "at least two")
})

test_that("compare_nicher rejects non-nicher objects", {
  expect_error(compare_nicher("a", "b"), "must be.*nicher")
})

test_that("compare_nicher rejects mismatched nobs", {
  fit_a <- .cv_cov_fit_po(n = 20L, seed = 1L)
  fit_b <- .cv_cov_fit_po(n = 15L, seed = 2L)
  expect_error(compare_nicher(fit_a, fit_b), "different.*nobs")
})

test_that("compare_nicher warns when mixing weighted and PO families", {
  fit_po <- .cv_cov_fit_po(seed = 1L)
  fit_w  <- .cv_cov_fit_w(seed = 1L)
  expect_warning(
    cmp <- compare_nicher(po = fit_po, w = fit_w),
    "Mixing fits"
  )
  expect_s3_class(cmp, "data.frame")
})

# ---- .env_fingerprint ----------------------------------------------------

test_that(".env_fingerprint returns NULL for NULL input", {
  expect_null(nicher:::.env_fingerprint(NULL))
})

test_that(".env_fingerprint produces stable fingerprints", {
  fp1 <- nicher:::.env_fingerprint(example_env_occ_2d)
  fp2 <- nicher:::.env_fingerprint(example_env_occ_2d)
  expect_identical(fp1, fp2)
  expect_equal(fp1$nrow, nrow(example_env_occ_2d))
  expect_equal(fp1$ncol, ncol(example_env_occ_2d))
})
