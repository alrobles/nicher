# Additional coverage tests for predict_nicher.R

skip_if_not_installed("terra")

test_that("predict.nicher rejects non-nicher input", {
  expect_error(predict.nicher("not_nicher", terra::rast()),
               "must be a.*nicher.*object")
})

test_that("predict.nicher rejects non-SpatRaster env", {
  set.seed(1)
  fit <- optimize_niche(
    env_occ = example_env_occ_2d, env_m = NULL,
    num_starts = 2L, likelihood = "presence_only",
    control = list(maxeval = 80L), verbose = FALSE
  )
  expect_error(predict(fit, as.matrix(example_env_occ_2d)),
               "must be a terra SpatRaster")
})

test_that("predict.nicher rejects malformed theta", {
  fit <- list(best = list(theta = c(1)), likelihood = "presence_only",
              var_names = "x")
  class(fit) <- "nicher"
  env <- terra::rast(nrows = 5, ncols = 5)
  terra::values(env) <- rnorm(25)
  names(env) <- "x"
  expect_error(predict(fit, env), "malformed|Cannot infer")
})

test_that("predict.nicher return_log gives log-suitability", {
  set.seed(2)
  fit <- optimize_niche(
    env_occ = example_env_occ_2d, env_m = NULL,
    num_starts = 2L, likelihood = "presence_only",
    control = list(maxeval = 80L), verbose = FALSE
  )
  r1 <- terra::rast(nrows = 10, ncols = 10)
  r2 <- terra::rast(nrows = 10, ncols = 10)
  terra::values(r1) <- rnorm(100, mean = 6, sd = 3)
  terra::values(r2) <- rnorm(100, mean = 200, sd = 100)
  env <- c(r1, r2)
  names(env) <- colnames(example_env_occ_2d)

  s_lin <- predict(fit, env, return_log = FALSE)
  s_log <- predict(fit, env, return_log = TRUE)
  v_lin <- terra::values(s_lin, na.rm = TRUE)
  v_log <- terra::values(s_log, na.rm = TRUE)

  expect_true(all(v_lin > 0 & v_lin <= 1))
  expect_true(all(v_log <= 0))
  expect_equal(as.numeric(log(v_lin)), as.numeric(v_log), tolerance = 1e-10)
})

test_that("predict.nicher handles unnamed layers with correct count", {
  set.seed(3)
  fit <- optimize_niche(
    env_occ = example_env_occ_2d, env_m = NULL,
    num_starts = 2L, likelihood = "presence_only",
    control = list(maxeval = 80L), verbose = FALSE
  )
  # Force var_names to NULL to test the positional fallback
  fit$var_names <- NULL
  r1 <- terra::rast(nrows = 5, ncols = 5, vals = rnorm(25, 6, 3))
  r2 <- terra::rast(nrows = 5, ncols = 5, vals = rnorm(25, 200, 100))
  env <- c(r1, r2)
  s <- predict(fit, env)
  expect_s4_class(s, "SpatRaster")
})

test_that(".reorder_env_for_predict errors on layer count mismatch", {
  r1 <- terra::rast(nrows = 5, ncols = 5, vals = rnorm(25))
  expect_error(
    nicher:::.reorder_env_for_predict(r1, var_names = NULL, p = 2L),
    "layers but the fitted model expects"
  )
})

test_that(".reorder_env_for_predict errors on unnamed layers when names expected", {
  r1 <- terra::rast(nrows = 5, ncols = 5, vals = rnorm(25))
  r2 <- terra::rast(nrows = 5, ncols = 5, vals = rnorm(25))
  env <- c(r1, r2)
  names(env) <- c("", "")
  expect_error(
    nicher:::.reorder_env_for_predict(env, var_names = c("bio1", "bio12"), p = 2L),
    "unnamed layers"
  )
})
