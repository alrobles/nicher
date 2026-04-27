# R/zzz.R
# Conditionally register the autoplot.nicher S3 method on ggplot2's
# `autoplot` generic. ggplot2 is a Suggests dependency; if it is not
# installed at load time, the method is simply not registered. This is
# the standard recipe (used by broom, ggdist, ggfortify, etc.) for
# wiring an S3 method into a Suggests-only generic without forcing
# every nicher install to drag in ggplot2.

.onLoad <- function(libname, pkgname) {
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    registerS3method(
      "autoplot", "nicher", autoplot.nicher,
      envir = asNamespace("ggplot2")
    )
  }
}
