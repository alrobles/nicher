test_that("nicher() returns an object of class 'nicher' for presence_only model", {
  fit <- nicher(spOccPnts, samMPts, model = "presence_only", itnmax = 20)
  expect_s3_class(fit, "nicher")
})

test_that("nicher object has required fields", {
  fit <- nicher(spOccPnts, samMPts, model = "presence_only", itnmax = 20)
  expect_true(all(c("loglik", "math_params", "bioscale_params",
                    "optimx_table", "model", "method") %in% names(fit)))
})

test_that("nicher model field is set correctly", {
  fit <- nicher(spOccPnts, samMPts, model = "presence_only", itnmax = 20)
  expect_equal(fit$model, "presence_only")
  expect_equal(fit$method, "mle")
})

test_that("nicher bioscale_params has mu, L, S, R, variances", {
  fit <- nicher(spOccPnts, samMPts, model = "presence_only", itnmax = 20)
  expect_true(all(c("mu", "L", "S", "R", "variances") %in% names(fit$bioscale_params)))
})

test_that("nicher bioscale_params: R is a correlation matrix", {
  fit <- nicher(spOccPnts, samMPts, model = "presence_only", itnmax = 20)
  R <- fit$bioscale_params$R
  expect_equal(diag(R), rep(1, nrow(R)), tolerance = 1e-8)
  expect_true(isSymmetric(R))
  expect_true(all(R >= -1 - 1e-8 & R <= 1 + 1e-8))
})

test_that("nicher bioscale_params: S = L %*% t(L)", {
  fit <- nicher(spOccPnts, samMPts, model = "presence_only", itnmax = 20)
  L <- fit$bioscale_params$L
  S <- fit$bioscale_params$S
  expect_equal(S, L %*% t(L), tolerance = 1e-8)
})

test_that("nicher loglik is a finite scalar", {
  fit <- nicher(spOccPnts, samMPts, model = "presence_only", itnmax = 20)
  expect_length(fit$loglik, 1)
  expect_true(is.finite(fit$loglik))
})

test_that("nicher() works for unweighted model", {
  fit <- nicher(spOccPnts, samMPts, model = "unweighted", itnmax = 20)
  expect_s3_class(fit, "nicher")
  expect_equal(fit$model, "unweighted")
})

test_that("nicher() works for weighted model", {
  fit <- nicher(spOccPnts, samMPts, model = "weighted", itnmax = 20)
  expect_s3_class(fit, "nicher")
  expect_equal(fit$model, "weighted")
})

test_that("print.nicher produces output without error", {
  fit <- nicher(spOccPnts, samMPts, model = "presence_only", itnmax = 20)
  expect_output(print(fit), "nicher")
  expect_output(print(fit), "Log-likelihood")
})

test_that("summary.nicher produces output and returns invisibly", {
  fit <- nicher(spOccPnts, samMPts, model = "presence_only", itnmax = 20)
  expect_output(s <- summary(fit), "nicher model summary")
  expect_s3_class(s, "summary.nicher")
  expect_true(all(c("loglik", "mu", "S", "R") %in% names(s)))
})

test_that("check_nicher raises error for non-nicher object", {
  expect_error(check_nicher(list()), "class 'nicher'")
})

test_that("check_nicher raises error for nicher with missing fields", {
  fit <- nicher(spOccPnts, samMPts, model = "presence_only", itnmax = 20)
  bad <- unclass(fit)
  bad$loglik <- NULL
  class(bad) <- "nicher"
  expect_error(check_nicher(bad), "missing required fields")
})
