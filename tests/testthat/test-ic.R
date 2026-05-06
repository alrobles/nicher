# Information criteria (AIC/BIC) and compare_nicher()

# Cheap reusable fixtures: IC helpers only need nicher object structure.
.mock_fit_ic <- function(likelihood = "weighted", loglik = -100, df = 5L,
                         env_occ = example_env_occ_2d,
                         env_m = example_env_m_2d) {
  theta <- stats::setNames(seq_len(df), paste0("theta", seq_len(df)))
  best <- list(
    theta = theta,
    loglik = loglik,
    loglik_unpenalised = loglik,
    convergence = 1L,
    nobs = nrow(env_occ),
    env_occ_fingerprint = nicher:::.env_fingerprint(env_occ),
    env_m_fingerprint = if (likelihood %in% c("presence_only", "skew_normal",
                                              "skew_t")) {
      NULL
    } else {
      nicher:::.env_fingerprint(env_m)
    }
  )
  solutions <- data.frame(
    start_id = 1L, loglik = loglik, convergence = 1L,
    stringsAsFactors = FALSE
  )
  nicher:::new_nicher(
    solutions = solutions,
    best = best,
    likelihood = likelihood,
    n_starts = 1L,
    var_names = colnames(env_occ)
  )
}
.fit_po <- function() .mock_fit_ic("presence_only", loglik = -120, df = 5L,
                                  env_m = NULL)
.fit_w <- function() .mock_fit_ic("weighted", loglik = -100, df = 5L)
.fit_skn_w <- function() .mock_fit_ic("skew_normal_weighted", loglik = -95,
                                      df = 7L)

# -----------------------------------------------------------------------
test_that("logLik.nicher returns un-penalised log-likelihood with df, nobs", {
  fit <- .fit_po()
  ll  <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_true(is.finite(as.numeric(ll)))
  expect_identical(attr(ll, "df"),   length(fit$best$theta))
  expect_identical(attr(ll, "nobs"), nrow(example_env_occ_2d))
  expect_equal(as.numeric(ll), fit$best$loglik, tolerance = 1e-8)
})

test_that("logLik.nicher equals -1 * neg_loglik wrapper at zero penalty", {
  fit <- .fit_w()
  ll  <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_identical(attr(ll, "df"),   length(fit$best$theta))
  expect_identical(attr(ll, "nobs"), nrow(example_env_occ_2d))
  expect_lte(fit$best$loglik, as.numeric(ll) + 1e-8)
})

# -----------------------------------------------------------------------
test_that("AIC.nicher and BIC.nicher follow exact formulas", {
  fit <- .fit_po()
  ll  <- as.numeric(logLik(fit))
  k   <- length(fit$best$theta)
  n   <- nobs(fit)
  expect_equal(AIC(fit), -2 * ll + 2 * k, tolerance = 1e-10)
  expect_equal(BIC(fit), -2 * ll + log(n) * k, tolerance = 1e-10)
  # Custom k recovers BIC
  expect_equal(AIC(fit, k = log(n)), BIC(fit), tolerance = 1e-10)
})

test_that("nobs.nicher returns nrow(env_occ)", {
  fit <- .fit_po()
  expect_identical(nobs(fit), nrow(example_env_occ_2d))
})

# -----------------------------------------------------------------------
test_that("compare_nicher produces a sensible IC table", {
  fit_w   <- .fit_w()
  fit_skn <- .fit_skn_w()
  cmp <- compare_nicher(weighted = fit_w, skew_normal_w = fit_skn)
  expect_s3_class(cmp, "data.frame")
  expect_setequal(
    names(cmp),
    c("model", "likelihood", "loglik", "df", "nobs",
      "AIC", "dAIC", "BIC", "dBIC", "weight_AIC", "convergence")
  )
  expect_equal(nrow(cmp), 2L)
  # Best AIC at top, dAIC[1] == 0
  expect_equal(cmp$dAIC[1L], 0)
  expect_true(all(cmp$dAIC >= -1e-8))
  expect_true(all(cmp$dBIC >= -1e-8))
  # Akaike weights sum to 1
  expect_equal(sum(cmp$weight_AIC), 1, tolerance = 1e-8)
  # Sample sizes match
  expect_true(all(cmp$nobs == nrow(example_env_occ_2d)))
})

test_that("compare_nicher refuses different env_occ", {
  fit1 <- .fit_w()
  occ2 <- example_env_occ_2d
  occ2[1L, 1L] <- occ2[1L, 1L] + 100
  fit2 <- .mock_fit_ic(env_occ = occ2)
  expect_error(
    compare_nicher(a = fit1, b = fit2),
    regexp = "different `env_occ`"
  )
})

test_that("compare_nicher refuses different env_m", {
  fit1 <- .fit_w()
  m2 <- example_env_m_2d
  m2[1L, 1L] <- m2[1L, 1L] + 100
  fit2 <- .mock_fit_ic(env_m = m2)
  expect_error(
    compare_nicher(a = fit1, b = fit2),
    regexp = "different `env_m`"
  )
})

test_that("compare_nicher warns when mixing PO with weighted", {
  fit_po <- .fit_po()
  fit_w  <- .fit_w()
  expect_warning(
    cmp <- compare_nicher(po = fit_po, weighted = fit_w),
    regexp = "Mixing fits with and without `env_m`"
  )
  expect_equal(nrow(cmp), 2L)
})

test_that("compare_nicher accepts a single named list", {
  fit_w   <- .fit_w()
  fit_skn <- .fit_skn_w()
  cmp_dots <- compare_nicher(weighted = fit_w, skn = fit_skn)
  cmp_list <- compare_nicher(list(weighted = fit_w, skn = fit_skn))
  # Same content modulo row order (sort_by default = "AIC" deterministic)
  expect_equal(cmp_dots$AIC, cmp_list$AIC)
  expect_equal(cmp_dots$model, cmp_list$model)
})

test_that("compare_nicher requires >= 2 nicher objects", {
  fit <- .fit_po()
  expect_error(compare_nicher(fit), regexp = "at least two")
  expect_error(
    compare_nicher(list(a = fit)),
    regexp = "at least two"
  )
})
