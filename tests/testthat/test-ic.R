# Information criteria (AIC/BIC) and compare_nicher()

# Cheap reusable fixtures: small num_starts, default breadth, fixed seed.
.fit_po <- function(seed = 1L) {
  set.seed(seed)
  optimize_niche(
    env_occ    = example_env_occ_2d,
    env_m      = NULL,
    num_starts = 5L,
    breadth    = 0.1,
    likelihood = "presence_only",
    verbose    = FALSE
  )
}
.fit_w <- function(seed = 1L) {
  set.seed(seed)
  optimize_niche(
    env_occ    = example_env_occ_2d,
    env_m      = example_env_m_2d,
    num_starts = 5L,
    breadth    = 0.1,
    likelihood = "weighted",
    verbose    = FALSE
  )
}
.fit_skn_w <- function(seed = 1L) {
  set.seed(seed)
  optimize_niche(
    env_occ    = example_env_occ_2d,
    env_m      = example_env_m_2d,
    num_starts = 5L,
    breadth    = 0.45,
    likelihood = "skew_normal_weighted",
    verbose    = FALSE
  )
}

# -----------------------------------------------------------------------
test_that("logLik.nicher returns un-penalised log-likelihood with df, nobs", {
  fit <- .fit_po()
  ll  <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_true(is.finite(as.numeric(ll)))
  expect_identical(attr(ll, "df"),   length(fit$best$theta))
  expect_identical(attr(ll, "nobs"), nrow(example_env_occ_2d))
  # presence_only has no penalty: best$loglik == loglik_unpenalised
  expect_equal(as.numeric(ll), fit$best$loglik, tolerance = 1e-8)
})

test_that("logLik.nicher equals -1 * neg_loglik wrapper at zero penalty", {
  fit <- .fit_w()
  ll  <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_identical(attr(ll, "df"),   length(fit$best$theta))
  expect_identical(attr(ll, "nobs"), nrow(example_env_occ_2d))
  # For penalised fits, best$loglik (penalised) >= loglik_unpenalised:
  # the penalty is non-negative, the optimiser MINIMISES (-loglik + pen),
  # so the stored best$loglik == -(neg_loglik + penalty) <= -neg_loglik
  # == loglik_unpenalised.  Equivalently:
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
  # Build a second fit with a perturbed env_occ -> different fingerprint.
  set.seed(2L)
  occ2 <- example_env_occ_2d
  occ2[1L, 1L] <- occ2[1L, 1L] + 100
  fit2 <- optimize_niche(
    env_occ = occ2, env_m = example_env_m_2d,
    num_starts = 5L, breadth = 0.1, likelihood = "weighted", verbose = FALSE
  )
  expect_error(
    compare_nicher(a = fit1, b = fit2),
    regexp = "different `env_occ`"
  )
})

test_that("compare_nicher refuses different env_m", {
  fit1 <- .fit_w()
  set.seed(3L)
  m2 <- example_env_m_2d
  m2[1L, 1L] <- m2[1L, 1L] + 100
  fit2 <- optimize_niche(
    env_occ = example_env_occ_2d, env_m = m2,
    num_starts = 5L, breadth = 0.1, likelihood = "weighted", verbose = FALSE
  )
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
