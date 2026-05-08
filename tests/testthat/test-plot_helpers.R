# Tests for plot_helpers.R internal functions.

test_that(".recover_mu_sigma returns correct geometry for Gaussian fit", {
  set.seed(1)
  fit <- optimize_niche(
    env_occ = example_env_occ_2d, env_m = NULL,
    num_starts = 2L, likelihood = "presence_only",
    control = list(maxeval = 80L), verbose = FALSE
  )
  pars <- nicher:::.recover_mu_sigma(fit)
  expect_equal(pars$p, 2L)
  expect_length(pars$mu, 2L)
  expect_equal(dim(pars$Sigma), c(2L, 2L))
  expect_true(all(eigen(pars$Sigma)$values > 0))
  expect_false(pars$is_skew)
  expect_false(pars$is_skew_t)
  expect_equal(pars$eta, 1.0)
})

test_that(".recover_mu_sigma works for skew_normal fit", {
  set.seed(2)
  fit <- optimize_niche(
    env_occ = example_env_occ_2d, env_m = NULL,
    num_starts = 2L, likelihood = "skew_normal",
    control = list(maxeval = 80L), verbose = FALSE
  )
  pars <- nicher:::.recover_mu_sigma(fit)
  expect_equal(pars$p, 2L)
  expect_true(pars$is_skew)
  expect_false(pars$is_skew_t)
  expect_length(pars$alpha, 2L)
})

test_that(".recover_mu_sigma errors on malformed theta", {
  fake <- list(best = list(theta = c(1)), likelihood = "presence_only")
  class(fake) <- "nicher"
  expect_error(nicher:::.recover_mu_sigma(fake), "missing or malformed|Cannot infer p")
})

test_that(".assert_nicher_2d accepts 2D fits and rejects non-nicher", {
  set.seed(3)
  fit <- optimize_niche(
    env_occ = example_env_occ_2d, env_m = NULL,
    num_starts = 2L, likelihood = "presence_only",
    control = list(maxeval = 80L), verbose = FALSE
  )
  expect_true(nicher:::.assert_nicher_2d(fit))
  expect_error(nicher:::.assert_nicher_2d("not_nicher"), "must be a nicher")
})

test_that(".assert_nicher_2d rejects 3D fits", {
  set.seed(4)
  fit3 <- optimize_niche(
    env_occ = example_env_occ_3d, env_m = NULL,
    num_starts = 2L, likelihood = "presence_only",
    control = list(maxeval = 80L), verbose = FALSE
  )
  expect_error(nicher:::.assert_nicher_2d(fit3), "must be a 2-D")
})

test_that(".coerce_xy returns a 2-col data.frame with correct names", {
  m <- matrix(1:6, ncol = 3)
  out <- nicher:::.coerce_xy(m, var_names = c("bio1", "bio12", "bio3"))
  expect_s3_class(out, "data.frame")
  expect_equal(ncol(out), 2L)
  expect_equal(names(out), c("bio1", "bio12"))

  out2 <- nicher:::.coerce_xy(m)
  expect_equal(names(out2), c("x1", "x2"))
})

test_that(".coerce_xy errors on < 2 columns", {
  expect_error(nicher:::.coerce_xy(matrix(1:3, ncol = 1)), "at least 2 columns")
  expect_error(nicher:::.coerce_xy("not_a_df"), "must be a matrix")
})

test_that(".ellipse_path returns a closed path for suitability level", {
  mu <- c(0, 0)
  Sigma <- diag(2)
  path <- nicher:::.ellipse_path(mu, Sigma, level = 0.5,
                                  level_type = "suitability", n = 50L)
  expect_s3_class(path, "data.frame")
  expect_equal(ncol(path), 2L)
  expect_equal(nrow(path), 51L)
  expect_equal(as.numeric(path[1, ]), as.numeric(path[51, ]))
})

test_that(".ellipse_path level = 1 gives a degenerate point", {
  path <- nicher:::.ellipse_path(c(0, 0), diag(2), level = 1,
                                  level_type = "suitability", n = 10L)
  expect_equal(nrow(path), 11L)
  expect_true(all(abs(path$x) < 1e-12))
  expect_true(all(abs(path$y) < 1e-12))
})

test_that(".ellipse_path supports chisq level_type", {
  path <- nicher:::.ellipse_path(c(1, 2), diag(2), level = 0.95,
                                  level_type = "chisq", n = 50L)
  expect_s3_class(path, "data.frame")
  expect_equal(nrow(path), 51L)
})

test_that(".ellipse_path warns on singular Sigma", {
  expect_warning(
    path <- nicher:::.ellipse_path(c(0, 0), matrix(0, 2, 2), level = 0.5),
    "rank-deficient"
  )
  expect_equal(nrow(path), 0L)
})

test_that(".ellipse_path rejects bad level values", {
  expect_error(nicher:::.ellipse_path(c(0, 0), diag(2), level = 0),
               "must be a single number")
  expect_error(nicher:::.ellipse_path(c(0, 0), diag(2), level = 1.5),
               "must be a single number")
})

test_that(".ellipse_paths stacks multiple levels", {
  mu <- c(0, 0)
  Sigma <- diag(2)
  paths <- nicher:::.ellipse_paths(mu, Sigma, level = c(0.5, 0.75),
                                    var_names = c("bio1", "bio12"), n = 20L)
  expect_s3_class(paths, "data.frame")
  expect_true("level" %in% names(paths))
  expect_equal(sort(unique(paths$level)), c(0.5, 0.75))
  expect_equal(names(paths)[1:2], c("bio1", "bio12"))
})

test_that(".ellipse_paths errors on bad level vector", {
  expect_error(
    nicher:::.ellipse_paths(c(0, 0), diag(2), level = c(-1, 0.5)),
    "values in \\(0, 1\\]"
  )
})

test_that(".assert_var_names_compatible checks model compatibility", {
  set.seed(5)
  fit_a <- optimize_niche(
    env_occ = example_env_occ_2d, env_m = NULL,
    num_starts = 2L, likelihood = "presence_only",
    control = list(maxeval = 80L), verbose = FALSE
  )
  fit_b <- optimize_niche(
    env_occ = example_env_occ_2d, env_m = NULL,
    num_starts = 2L, likelihood = "presence_only",
    control = list(maxeval = 80L), verbose = FALSE
  )
  expect_invisible(
    nicher:::.assert_var_names_compatible(list(fit_a, fit_b))
  )
  expect_error(
    nicher:::.assert_var_names_compatible(list()),
    "non-empty list"
  )
})

test_that(".linewidth_param returns the right parameter name", {
  lw <- nicher:::.linewidth_param(1.5)
  expect_true("linewidth" %in% names(lw) || "size" %in% names(lw))
  val <- if (!is.null(lw$linewidth)) lw$linewidth else lw$size
  expect_equal(val, 1.5)
})

test_that(".require_ggplot2 errors when ggplot2 is absent", {
  skip_if(requireNamespace("ggplot2", quietly = TRUE),
          message = "ggplot2 installed; cannot test absence path")
  expect_error(nicher:::.require_ggplot2("test_fn"), "ggplot2")
})

test_that(".gl_quad_20 returns 20 nodes and weights", {
  gl <- nicher:::.gl_quad_20()
  expect_length(gl$nodes, 20L)
  expect_length(gl$weights, 20L)
  expect_true(all(gl$nodes > 0))
  expect_true(all(gl$weights > 0))
})
