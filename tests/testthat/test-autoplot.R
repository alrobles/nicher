# Tests for autoplot.nicher() and the composable layer constructors.

skip_if_no_ggplot2 <- function() {
  testthat::skip_if_not_installed("ggplot2")
}

fit_2d <- function(num_starts = 5L) {
  set.seed(2024L)
  optimize_niche(
    env_occ    = example_env_occ_2d,
    env_m      = example_env_m_2d,
    num_starts = num_starts
  )
}

fit_3d <- function(num_starts = 3L) {
  set.seed(2024L)
  optimize_niche(
    env_occ    = example_env_occ_3d,
    env_m      = example_env_m_3d,
    num_starts = num_starts
  )
}

test_that("autoplot.nicher returns a ggplot object with 3 layers", {
  skip_if_no_ggplot2()
  fit <- fit_2d()
  p <- ggplot2::autoplot(fit,
                          env_occ = example_env_occ_2d,
                          env_m   = example_env_m_2d)
  expect_s3_class(p, "ggplot")
  expect_s3_class(p, "gg")
  expect_length(p$layers, 3L)
  expect_identical(p$labels$x, fit$var_names[1])
  expect_identical(p$labels$y, fit$var_names[2])
})

test_that("autoplot.nicher rejects 3-D fits with an informative error", {
  skip_if_no_ggplot2()
  fit3 <- fit_3d()
  expect_error(
    ggplot2::autoplot(fit3,
                       env_occ = example_env_occ_3d,
                       env_m   = example_env_m_3d),
    "must be a 2-D nicher fit"
  )
})

test_that("autoplot.nicher requires explicit env_occ and env_m", {
  skip_if_no_ggplot2()
  fit <- fit_2d()
  expect_error(ggplot2::autoplot(fit), "must be provided explicitly")
  expect_error(
    ggplot2::autoplot(fit, env_occ = example_env_occ_2d),
    "must be provided explicitly"
  )
  expect_error(
    ggplot2::autoplot(fit, env_m = example_env_m_2d),
    "must be provided explicitly"
  )
})

test_that("layer constructors are composable on a fresh ggplot", {
  skip_if_no_ggplot2()
  fit <- fit_2d()
  p <- ggplot2::ggplot() +
    geom_nicher_background(example_env_m_2d) +
    geom_nicher_occ(example_env_occ_2d) +
    geom_nicher_ellipse(fit)
  expect_s3_class(p, "ggplot")
  expect_length(p$layers, 3L)
  # Building the plot should not error.
  expect_silent(suppressMessages(ggplot2::ggplot_build(p)))
})

test_that("geom_nicher_ellipse with multiple levels yields a single layer", {
  skip_if_no_ggplot2()
  fit <- fit_2d()
  layer <- geom_nicher_ellipse(fit, level = c(0.95, 0.5))
  expect_s3_class(layer, "Layer")
  expect_true(nrow(layer$data) > 0L)
  # 2 closed ellipse paths => 2 unique levels in the layer's data.
  expect_setequal(unique(layer$data$level), c(0.95, 0.5))
})

test_that("matrix and data.frame env_occ/env_m both work", {
  skip_if_no_ggplot2()
  fit <- fit_2d()
  m_mat <- as.matrix(example_env_m_2d)
  o_mat <- as.matrix(example_env_occ_2d)
  expect_s3_class(
    ggplot2::autoplot(fit, env_occ = o_mat, env_m = m_mat),
    "ggplot"
  )
  expect_s3_class(
    ggplot2::autoplot(fit, env_occ = as.data.frame(o_mat),
                       env_m = as.data.frame(m_mat)),
    "ggplot"
  )
})

test_that("iso-suitability geometry: every point on level=alpha satisfies S(x) = alpha", {
  skip_if_no_ggplot2()
  fit <- fit_2d()
  ms <- nicher:::.recover_mu_sigma(fit)
  for (lv in c(0.05, 0.5, 0.95)) {
    pth <- nicher:::.ellipse_path(ms$mu, ms$Sigma, lv,
                                   level_type = "suitability",
                                   n = 50L)
    S <- apply(as.matrix(pth), 1L, function(x) {
      d2 <- as.numeric(t(x - ms$mu) %*% solve(ms$Sigma) %*% (x - ms$mu))
      exp(-0.5 * d2)
    })
    expect_true(max(abs(S - lv)) < 1e-9,
                info = sprintf("level = %g", lv))
  }
})

test_that("chisq geometry: every point on level=q satisfies d^2 = qchisq(q, 2)", {
  skip_if_no_ggplot2()
  fit <- fit_2d()
  ms <- nicher:::.recover_mu_sigma(fit)
  for (q in c(0.5, 0.9, 0.95)) {
    pth <- nicher:::.ellipse_path(ms$mu, ms$Sigma, q,
                                   level_type = "chisq",
                                   n = 50L)
    d2 <- apply(as.matrix(pth), 1L, function(x) {
      as.numeric(t(x - ms$mu) %*% solve(ms$Sigma) %*% (x - ms$mu))
    })
    expect_true(max(abs(d2 - stats::qchisq(q, df = 2L))) < 1e-9,
                info = sprintf("q = %g", q))
  }
})

test_that("nicher_compare_plot accepts shared env_occ and returns a ggplot", {
  skip_if_no_ggplot2()
  fit_a <- fit_2d(num_starts = 4L)
  fit_b <- fit_2d(num_starts = 6L)
  p <- nicher_compare_plot(
    models  = list(a = fit_a, b = fit_b),
    env_occ = example_env_occ_2d,
    env_m   = example_env_m_2d
  )
  expect_s3_class(p, "ggplot")
  # background + occurrences (shared) + ellipses (combined) = 3
  expect_length(p$layers, 3L)
})

test_that("nicher_compare_plot accepts a named-list env_occ", {
  skip_if_no_ggplot2()
  fit_a <- fit_2d(num_starts = 4L)
  fit_b <- fit_2d(num_starts = 6L)
  half <- floor(nrow(example_env_occ_2d) / 2L)
  occ_a <- example_env_occ_2d[seq_len(half), , drop = FALSE]
  occ_b <- example_env_occ_2d[(half + 1L):nrow(example_env_occ_2d), ,
                               drop = FALSE]
  p <- nicher_compare_plot(
    models  = list(a = fit_a, b = fit_b),
    env_occ = list(a = occ_a, b = occ_b),
    env_m   = example_env_m_2d
  )
  expect_s3_class(p, "ggplot")
})

test_that("nicher_compare_plot rejects unnamed models", {
  skip_if_no_ggplot2()
  fit_a <- fit_2d(num_starts = 4L)
  fit_b <- fit_2d(num_starts = 4L)
  expect_error(
    nicher_compare_plot(
      models  = list(fit_a, fit_b),
      env_occ = example_env_occ_2d,
      env_m   = example_env_m_2d
    ),
    "named list"
  )
})

test_that("nicher_compare_plot rejects var_names mismatch across models", {
  skip_if_no_ggplot2()
  fit_a <- fit_2d(num_starts = 4L)
  fit_b <- fit_2d(num_starts = 4L)
  fit_b$var_names <- c("renamed_x", "renamed_y")
  expect_error(
    nicher_compare_plot(
      models  = list(a = fit_a, b = fit_b),
      env_occ = example_env_occ_2d,
      env_m   = example_env_m_2d
    ),
    "Models are not compatible"
  )
})

test_that("nicher_compare_plot rejects mismatched env_occ list names", {
  skip_if_no_ggplot2()
  fit_a <- fit_2d(num_starts = 4L)
  fit_b <- fit_2d(num_starts = 4L)
  expect_error(
    nicher_compare_plot(
      models  = list(a = fit_a, b = fit_b),
      env_occ = list(a = example_env_occ_2d, c = example_env_occ_2d),
      env_m   = example_env_m_2d
    ),
    "match `names\\(models\\)`"
  )
})

test_that("layers with too-few env columns error", {
  skip_if_no_ggplot2()
  expect_error(
    geom_nicher_background(matrix(1:5, ncol = 1)),
    "must have at least 2 columns"
  )
  expect_error(
    geom_nicher_occ(matrix(1:5, ncol = 1)),
    "must have at least 2 columns"
  )
})

test_that("legacy nicher fit without var_names still plots", {
  skip_if_no_ggplot2()
  fit <- fit_2d()
  fit$var_names <- NULL
  p <- ggplot2::autoplot(fit,
                          env_occ = example_env_occ_2d,
                          env_m   = example_env_m_2d)
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$x, "x1")
  expect_identical(p$labels$y, "x2")
})
