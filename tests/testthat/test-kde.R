test_that("kde_bandwidth_silverman returns a symmetric positive-definite matrix", {
  H <- kde_bandwidth_silverman(as.matrix(samMPts))

  expect_true(is.matrix(H))
  expect_equal(nrow(H), ncol(H))
  expect_equal(nrow(H), ncol(samMPts))
  # Symmetric
  expect_equal(H, t(H))
  # Positive definite: all eigenvalues > 0
  expect_true(all(eigen(H, only.values = TRUE)$values > 0))
})

test_that("kde_bandwidth_silverman scales with sample size", {
  set.seed(42)
  small <- matrix(rnorm(20 * 2), ncol = 2)
  large <- matrix(rnorm(200 * 2), ncol = 2)
  H_small <- kde_bandwidth_silverman(small)
  H_large <- kde_bandwidth_silverman(large)
  # Larger sample => smaller bandwidth (Silverman rule h ~ n^{-2/(p+4)})
  expect_true(H_small[1, 1] > H_large[1, 1])
})

test_that("kde_eval_cpp returns a finite numeric vector of length m", {
  H   <- kde_bandwidth_silverman(as.matrix(samMPts))
  lkd <- kde_eval_cpp(as.matrix(spOccPnts), as.matrix(samMPts), H)

  expect_length(lkd, nrow(spOccPnts))
  expect_type(lkd, "double")
  expect_true(all(is.finite(lkd)))
})

test_that("kde_bandwidth_cached returns same result as kde_bandwidth_silverman", {
  kde_cache_clear()
  H_direct <- kde_bandwidth_silverman(as.matrix(samMPts))
  H_cached <- kde_bandwidth_cached(samMPts)

  expect_equal(H_cached, H_direct)
})

test_that("kde_bandwidth_cached reuses cached bandwidth on second call", {
  kde_cache_clear()
  H1 <- kde_bandwidth_cached(samMPts)
  H2 <- kde_bandwidth_cached(samMPts)

  expect_identical(H1, H2)
})

test_that("kde_cache_clear removes cached entries", {
  kde_bandwidth_cached(samMPts)  # populate cache
  kde_cache_clear()

  expect_equal(length(ls(nicher:::.kde_cache)), 0L)
})

test_that("kde_eval_cached matches kde_eval_cpp with explicit bandwidth", {
  H    <- kde_bandwidth_silverman(as.matrix(samMPts))
  lkd1 <- kde_eval_cpp(as.matrix(spOccPnts), as.matrix(samMPts), H)
  lkd2 <- kde_eval_cached(spOccPnts, samMPts, H = H)

  expect_equal(lkd1, lkd2)
})

test_that("kde_eval_cached auto-computes bandwidth when H is NULL", {
  kde_cache_clear()
  lkd <- kde_eval_cached(spOccPnts, samMPts)

  expect_length(lkd, nrow(spOccPnts))
  expect_true(all(is.finite(lkd)))
})
