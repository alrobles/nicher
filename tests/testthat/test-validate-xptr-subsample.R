test_that(".validate_xptr_result evaluates the same subsampled objective the optimizer used", {
  # Build an env_m large enough that the default cap (10000) forces
  # subsampling. If .validate_xptr_result evaluated on the FULL env_m it
  # would fire a spurious pointer-safety warning, since the optimizer ran
  # against the subsampled likelihood.
  set.seed(0)
  p     <- 2L
  occ   <- as.data.frame(matrix(stats::rnorm(80 * p), ncol = p,
                                dimnames = list(NULL, c("x1", "x2"))))
  M_big <- as.data.frame(matrix(stats::rnorm(12000 * p), ncol = p,
                                dimnames = list(NULL, c("x1", "x2"))))

  expect_silent(
    res <- optimize_niche(
      env_occ    = occ,
      env_m      = M_big,
      num_starts = 2L,
      breadth    = 0.1,
      likelihood = "weighted",
      backend    = "cpp",
      seed       = 1L,
      control    = list(maxeval = 200L)
    )
  )
  expect_s3_class(res, "nicher")
})
