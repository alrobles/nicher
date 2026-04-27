test_that(".resolve_weighted_inputs caps subsamples at 10000 by default", {
  set.seed(1)
  p <- 2L
  fake_m <- matrix(stats::rnorm(50000 * p), ncol = p)
  fake_occ <- fake_m[1:50, , drop = FALSE]
  colnames(fake_m)   <- c("x1", "x2")
  colnames(fake_occ) <- c("x1", "x2")

  res <- nicher:::.resolve_weighted_inputs(
    env_occ = as.data.frame(fake_occ),
    env_m   = as.data.frame(fake_m),
    seed    = 1L
  )
  expect_equal(length(res$den_idx), 10000L)
  expect_equal(length(res$kde_idx), 10000L)
  expect_length(res$w_occ, nrow(fake_occ))
  expect_length(res$w_den, 10000L)
})

test_that(".resolve_weighted_inputs uses all rows when n_m < cap", {
  set.seed(2)
  p <- 2L
  fake_m <- matrix(stats::rnorm(2000 * p), ncol = p)
  fake_occ <- fake_m[1:50, , drop = FALSE]
  colnames(fake_m)   <- c("x1", "x2")
  colnames(fake_occ) <- c("x1", "x2")

  res <- nicher:::.resolve_weighted_inputs(
    env_occ = as.data.frame(fake_occ),
    env_m   = as.data.frame(fake_m),
    seed    = 2L
  )
  expect_equal(length(res$den_idx), 2000L)
  expect_equal(length(res$kde_idx), 2000L)
})

test_that(".resolve_weighted_inputs warns below the representative floor", {
  set.seed(3)
  p <- 3L
  # floor = max(500, 50 * 2^3) = max(500, 400) = 500. Stay below.
  fake_m <- matrix(stats::rnorm(200 * p), ncol = p)
  fake_occ <- fake_m[1:20, , drop = FALSE]
  colnames(fake_m)   <- c("x1", "x2", "x3")
  colnames(fake_occ) <- c("x1", "x2", "x3")

  expect_warning(
    nicher:::.resolve_weighted_inputs(
      env_occ = as.data.frame(fake_occ),
      env_m   = as.data.frame(fake_m),
      seed    = 3L
    ),
    regexp = "below the recommended minimum"
  )
})

test_that(".resolve_weighted_inputs is deterministic given a seed", {
  set.seed(4)
  p <- 2L
  fake_m <- matrix(stats::rnorm(20000 * p), ncol = p)
  fake_occ <- fake_m[1:30, , drop = FALSE]
  colnames(fake_m)   <- c("x1", "x2")
  colnames(fake_occ) <- c("x1", "x2")

  r1 <- nicher:::.resolve_weighted_inputs(
    env_occ = as.data.frame(fake_occ),
    env_m   = as.data.frame(fake_m),
    seed    = 7L
  )
  r2 <- nicher:::.resolve_weighted_inputs(
    env_occ = as.data.frame(fake_occ),
    env_m   = as.data.frame(fake_m),
    seed    = 7L
  )
  expect_identical(r1$den_idx, r2$den_idx)
  expect_identical(r1$kde_idx, r2$kde_idx)
})
