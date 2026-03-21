# tests/testthat/test-niche-wrappers.R
#
# Tests for the single-start niche optimization wrappers:
#   niche_unweighted(), niche_presence_only(), niche_weighted()
#
# These tests use the built-in 3-D example datasets so that they run quickly
# on CRAN. Heavier multi-start scenarios are guarded with skip_on_cran().

# ---------------------------------------------------------------------------
# Shared test data
# ---------------------------------------------------------------------------

occ3 <- as.matrix(example_env_occ_3d)
M3   <- as.matrix(example_env_m_3d)

set.seed(42)
theta0_3d <- start_theta(example_env_occ_3d)

# Precomputed KDE weights for niche_weighted tests
set.seed(123)
den_idx3 <- sample.int(nrow(M3), 300L)
kde_idx3 <- sample.int(nrow(M3), 600L)
pre3     <- kde_gaussian(M3[den_idx3, ], M3[kde_idx3, ])

# ---------------------------------------------------------------------------
# Test 1 – niche_unweighted returns correct structure
# ---------------------------------------------------------------------------
test_that("niche_unweighted returns a list with value, conv, and theta", {
  res <- niche_unweighted(occ = occ3, M = M3, start = theta0_3d)

  expect_type(res, "list")
  expect_named(res, c("value", "conv", "theta"), ignore.order = TRUE)

  expect_true(is.numeric(res$value) && length(res$value) == 1L,
              label = "value is a scalar numeric")
  expect_true(is.integer(res$conv) && length(res$conv) == 1L,
              label = "conv is a scalar integer")
  expect_true(is.numeric(res$theta),
              label = "theta is a numeric vector")
  expect_equal(length(res$theta), length(theta0_3d),
               label = "theta has the same length as start")
})

# ---------------------------------------------------------------------------
# Test 2 – niche_unweighted returns a finite log-likelihood
# ---------------------------------------------------------------------------
test_that("niche_unweighted returns a finite log-likelihood", {
  res <- niche_unweighted(occ = occ3, M = M3, start = theta0_3d)

  expect_true(is.finite(res$value),
              label = "niche_unweighted value is finite")
})

# ---------------------------------------------------------------------------
# Test 3 – niche_presence_only returns correct structure
# ---------------------------------------------------------------------------
test_that("niche_presence_only returns a list with value, conv, and theta", {
  res <- niche_presence_only(occ = occ3, start = theta0_3d)

  expect_type(res, "list")
  expect_named(res, c("value", "conv", "theta"), ignore.order = TRUE)

  expect_true(is.numeric(res$value) && length(res$value) == 1L,
              label = "value is a scalar numeric")
  expect_true(is.integer(res$conv) && length(res$conv) == 1L,
              label = "conv is a scalar integer")
  expect_true(is.numeric(res$theta),
              label = "theta is a numeric vector")
  expect_equal(length(res$theta), length(theta0_3d),
               label = "theta has the same length as start")
})

# ---------------------------------------------------------------------------
# Test 4 – niche_presence_only returns a finite log-likelihood
# ---------------------------------------------------------------------------
test_that("niche_presence_only returns a finite log-likelihood", {
  res <- niche_presence_only(occ = occ3, start = theta0_3d)

  expect_true(is.finite(res$value),
              label = "niche_presence_only value is finite")
})

# ---------------------------------------------------------------------------
# Test 5 – niche_weighted returns correct structure
# ---------------------------------------------------------------------------
test_that("niche_weighted returns a list with value, conv, and theta", {
  res <- niche_weighted(
    occ           = occ3,
    M             = M3,
    den_idx       = den_idx3,
    kde_idx       = kde_idx3,
    precomp_w_den = pre3,
    start         = theta0_3d
  )

  expect_type(res, "list")
  expect_named(res, c("value", "conv", "theta"), ignore.order = TRUE)

  expect_true(is.numeric(res$value) && length(res$value) == 1L,
              label = "value is a scalar numeric")
  expect_true(is.integer(res$conv) && length(res$conv) == 1L,
              label = "conv is a scalar integer")
  expect_true(is.numeric(res$theta),
              label = "theta is a numeric vector")
  expect_equal(length(res$theta), length(theta0_3d),
               label = "theta has the same length as start")
})

# ---------------------------------------------------------------------------
# Test 6 – niche_weighted returns a finite log-likelihood
# ---------------------------------------------------------------------------
test_that("niche_weighted returns a finite log-likelihood", {
  res <- niche_weighted(
    occ           = occ3,
    M             = M3,
    den_idx       = den_idx3,
    kde_idx       = kde_idx3,
    precomp_w_den = pre3,
    start         = theta0_3d
  )

  expect_true(is.finite(res$value),
              label = "niche_weighted value is finite")
})

# ---------------------------------------------------------------------------
# Test 7 – precomp_w_den length mismatch raises an informative error
# ---------------------------------------------------------------------------
test_that("niche_weighted errors when precomp_w_den length does not match den_idx", {
  expect_error(
    niche_weighted(
      occ           = occ3,
      M             = M3,
      den_idx       = den_idx3,
      kde_idx       = kde_idx3,
      precomp_w_den = pre3[-1L],   # one element short
      start         = theta0_3d
    ),
    regexp = "precomp_w_den.*must have the same length",
    label  = "length mismatch produces informative error"
  )
})

# ---------------------------------------------------------------------------
# Test 8 – Non-finite starting values raise an error
# ---------------------------------------------------------------------------
test_that("wrappers error on non-finite starting values", {
  bad_start <- theta0_3d
  bad_start[1L] <- Inf

  expect_error(niche_unweighted(occ3, M3, bad_start),
               regexp = "finite",
               label  = "niche_unweighted errors on Inf start")

  expect_error(niche_presence_only(occ3, bad_start),
               regexp = "finite",
               label  = "niche_presence_only errors on Inf start")

  expect_error(
    niche_weighted(occ3, M3, den_idx3, kde_idx3, pre3, bad_start),
    regexp = "finite",
    label  = "niche_weighted errors on Inf start"
  )
})

# ---------------------------------------------------------------------------
# Test 9 – Multi-start loop: all three models run for every starting point
#           (off-CRAN only due to heavier computation with full Sobol starts)
# ---------------------------------------------------------------------------
test_that("multi-start loop completes without error for UN, PO, and W (3D)", {
  skip_on_cran()

  set.seed(2026)
  starts_df <- start_theta_multiple(
    env_data   = as.data.frame(occ3),
    num_starts = 5L,
    method     = "uniform"
  )
  starts <- lapply(seq_len(nrow(starts_df)), function(i) as.numeric(starts_df[i, ]))

  for (i in seq_along(starts)) {
    s <- starts[[i]]

    # Unweighted
    un <- niche_unweighted(occ = occ3, M = M3, start = s)
    expect_true(is.finite(un$value),
                label = sprintf("UN start %d: finite value", i))

    # Presence-only
    po <- niche_presence_only(occ = occ3, start = s)
    expect_true(is.finite(po$value),
                label = sprintf("PO start %d: finite value", i))

    # Weighted
    w <- niche_weighted(
      occ           = occ3,
      M             = M3,
      den_idx       = den_idx3,
      kde_idx       = kde_idx3,
      precomp_w_den = pre3,
      start         = s
    )
    expect_true(is.finite(w$value),
                label = sprintf("W start %d: finite value", i))
  }
})
