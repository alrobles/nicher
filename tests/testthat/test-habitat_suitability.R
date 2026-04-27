skip_if_not_installed("terra")

build_test_env <- function(seed = 1L, nr = 25L, nc = 30L) {
  set.seed(seed)
  r1 <- terra::rast(nrows = nr, ncols = nc, vals = stats::rnorm(nr * nc))
  r2 <- terra::rast(nrows = nr, ncols = nc, vals = stats::rnorm(nr * nc))
  out <- c(r1, r2)
  names(out) <- c("bio1", "bio12")
  out
}

test_that("habitat_suitability returns a one-layer SpatRaster in (0, 1]", {
  env <- build_test_env()
  param <- list(mu = c(0, 0), Sigma = diag(2L))
  s <- habitat_suitability(param, env)
  expect_s4_class(s, "SpatRaster")
  expect_equal(terra::nlyr(s), 1L)
  v <- terra::values(s, na.rm = TRUE)
  expect_true(all(v > 0 & v <= 1))
  expect_true("suitability" %in% names(s))
})

test_that("block-loop output matches the unblocked C++ kernel exactly", {
  env <- build_test_env(seed = 2L)
  mu <- c(0.2, -0.3)
  A <- matrix(c(1.5, 0.5, 0.5, 0.8), 2, 2); Sigma <- crossprod(A)
  s <- habitat_suitability(list(mu = mu, Sigma = Sigma), env)
  ref_mat <- as.matrix(terra::values(env))
  L_inv <- backsolve(t(chol(Sigma)), diag(2L), upper.tri = FALSE)
  ref <- niche_suitability_cpp(
    as.numeric(ref_mat), c(nrow(ref_mat), 2L), mu, L_inv, FALSE, 0L
  )
  got <- as.numeric(terra::values(s))
  expect_equal(got, ref, tolerance = 1e-14)
})

test_that("NA cells in env produce NA cells in the suitability map", {
  env <- build_test_env(seed = 3L)
  v1 <- terra::values(env[[1L]])
  v1[c(1L, 5L, 100L)] <- NA_real_
  terra::values(env[[1L]]) <- v1
  s <- habitat_suitability(list(mu = c(0, 0), Sigma = diag(2L)), env)
  v <- terra::values(s)
  expect_true(all(is.na(v[c(1L, 5L, 100L)])))
  expect_true(all(is.finite(v[-c(1L, 5L, 100L)])))
})

test_that("return_log = TRUE emits log-suitability layer name and <= 0 values", {
  env <- build_test_env(seed = 4L)
  s <- habitat_suitability(list(mu = c(0, 0), Sigma = diag(2L)),
                           env, return_log = TRUE)
  expect_true("log_suitability" %in% names(s))
  v <- terra::values(s, na.rm = TRUE)
  expect_true(all(v <= 0))
})

test_that("output writes a readable GeoTIFF with correct geometry", {
  env <- build_test_env(seed = 5L)
  tmp <- tempfile(fileext = ".tif")
  on.exit(unlink(tmp), add = TRUE)
  invisible(habitat_suitability(
    list(mu = c(0, 0), Sigma = diag(2L)),
    env,
    output    = tmp,
    overwrite = TRUE
  ))
  expect_true(file.exists(tmp))
  s <- terra::rast(tmp)
  expect_equal(terra::nlyr(s), 1L)
  expect_true(terra::compareGeom(s, env, stopOnError = FALSE))
})

test_that("degenerate Sigma yields a warning + NA raster, not an abort", {
  env <- build_test_env(seed = 6L)
  Sigma_singular <- matrix(c(1, 1, 1, 1), 2, 2)
  expect_warning(
    s <- habitat_suitability(
      list(mu = c(0, 0), Sigma = Sigma_singular),
      env
    ),
    "chol\\(Sigma\\) failed"
  )
  expect_true(all(is.na(terra::values(s))))
})

test_that("layer-count mismatch and bad inputs error clearly", {
  env <- build_test_env(seed = 7L)
  expect_error(
    habitat_suitability(list(mu = c(0, 0, 0), Sigma = diag(3L)), env),
    "one layer per environmental variable"
  )
  expect_error(habitat_suitability(list(mu = c(0, 0)), env),
               "subset")  # checkmate: missing Sigma
  expect_error(habitat_suitability(list(mu = c(0, 0), Sigma = diag(2L)),
                                    env = matrix(0, 4, 2)),
               "SpatRaster")
})
