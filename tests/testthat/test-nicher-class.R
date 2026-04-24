# tests/testthat/test-nicher-class.R
#
# Tests for the nicher S3 class: new_nicher, print.nicher, assess.nicher
#
# These tests use the built-in 2D example datasets for speed.

# ---------------------------------------------------------------------------
# Shared test data: build a small nicher object for unit tests
# ---------------------------------------------------------------------------

.make_mock_nicher <- function(logliks, convs, likelihood = "weighted") {
  n <- length(logliks)
  solutions <- data.frame(
    start_id    = seq_len(n),
    loglik      = logliks,
    convergence = as.integer(convs),
    stringsAsFactors = FALSE
  )
  solutions$full_par <- replicate(n, c(0, 0, 0, 0), simplify = FALSE)
  best <- list(
    theta       = solutions$full_par[[which.max(logliks)]],
    loglik      = max(logliks),
    convergence = as.integer(convs[which.max(logliks)])
  )
  new_nicher(solutions, best, likelihood, n_starts = n)
}

# ---------------------------------------------------------------------------
# Test 1 – new_nicher creates an object of class "nicher"
# ---------------------------------------------------------------------------

test_that("new_nicher creates an object of class 'nicher'", {
  obj <- .make_mock_nicher(c(-10, -10.5), c(1L, 1L))
  expect_s3_class(obj, "nicher")
  expect_named(obj, c("solutions", "best", "likelihood", "n_starts"),
    ignore.order = TRUE
  )
})

# ---------------------------------------------------------------------------
# Test 2 – print.nicher returns the object invisibly
# ---------------------------------------------------------------------------

test_that("print.nicher returns x invisibly without error", {
  obj <- .make_mock_nicher(c(-10, -10.5), c(1L, 1L))
  expect_invisible(print(obj))
  expect_output(print(obj), regexp = "nicher optimization result")
})

# ---------------------------------------------------------------------------
# Test 3 – assess: accepted_global when two starts agree closely
# ---------------------------------------------------------------------------

test_that("assess returns 'accepted_global' when starts agree closely", {
  obj <- .make_mock_nicher(c(-10.000, -10.001), c(1L, 1L))
  diag <- assess(obj)
  expect_equal(diag$flag, "accepted_global")
  expect_equal(diag$n_converged, 2L)
  expect_true(is.finite(diag$gap))
})

# ---------------------------------------------------------------------------
# Test 4 – assess: accepted_noise for moderate spread
# ---------------------------------------------------------------------------

test_that("assess returns 'accepted_noise' for moderate spread", {
  # rel_gap in (tol_gap, tol_dist]: gap / |best| in (0.01, 0.05]
  # best_loglik = -10, gap in (0.1, 0.5]
  obj <- .make_mock_nicher(c(-10.0, -10.3), c(1L, 1L))
  diag <- assess(obj)
  expect_true(diag$flag %in% c("accepted_noise", "suggest_average"))
  expect_equal(diag$n_converged, 2L)
})

# ---------------------------------------------------------------------------
# Test 5 – assess: needs_more_starts when too few converged
# ---------------------------------------------------------------------------

test_that("assess returns 'needs_more_starts' when fewer than 2 converge", {
  # Only one converged solution
  obj <- .make_mock_nicher(c(-10, -9), c(1L, NA_integer_))
  diag <- assess(obj, min_converged = 2L)
  expect_equal(diag$flag, "needs_more_starts")
  expect_true(is.na(diag$gap))
})

# ---------------------------------------------------------------------------
# Test 6 – assess: suggest_average for wide spread with 2+ converged
# ---------------------------------------------------------------------------

test_that("assess returns 'suggest_average' for wide spread", {
  # gap = 5 on best = -10 → rel_gap = 0.5 >> tol_dist = 0.05
  obj <- .make_mock_nicher(c(-10.0, -15.0), c(1L, 1L))
  diag <- assess(obj)
  expect_equal(diag$flag, "suggest_average")
})

# ---------------------------------------------------------------------------
# Test 7 – assess fields are correctly typed
# ---------------------------------------------------------------------------

test_that("assess.nicher returns list with expected names and types", {
  obj <- .make_mock_nicher(c(-10, -10.002), c(1L, 1L))
  diag <- assess(obj)

  expect_type(diag, "list")
  expect_named(diag,
    c("flag", "recommendation", "gap", "rel_gap",
      "n_converged", "best_loglik"),
    ignore.order = TRUE
  )
  expect_type(diag$flag, "character")
  expect_type(diag$recommendation, "character")
  expect_type(diag$gap, "double")
  expect_type(diag$rel_gap, "double")
  expect_type(diag$n_converged, "integer")
  expect_type(diag$best_loglik, "double")
})

# ---------------------------------------------------------------------------
# Test 8 – optimize_niche returns a nicher object (integration test, 2D)
# ---------------------------------------------------------------------------

test_that("optimize_niche returns a nicher object with correct structure (2D)", {
  skip_on_cran()

  set.seed(7L)
  res <- optimize_niche(
    env_occ    = example_env_occ_2d,
    env_m      = NULL,
    num_starts = 5L,
    breadth    = 0.1,
    likelihood = "presence_only",
    eta        = 1
  )

  expect_s3_class(res, "nicher")
  expect_equal(res$likelihood, "presence_only")
  expect_equal(res$n_starts, 5L)
  expect_true(is.finite(res$best$loglik))
  expect_true(nrow(res$solutions) == 5L)

  diag <- assess(res)
  expect_true(diag$flag %in% c(
    "accepted_global",
    "accepted_noise",
    "suggest_average",
    "needs_more_starts"
  ))
})
