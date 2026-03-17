# Internal KDE bandwidth cache: an environment keyed by a lightweight
# checksum of the reference sample so that the expensive Silverman
# bandwidth computation is performed only once per unique background sample.
.kde_cache <- new.env(parent = emptyenv())

#' Compute and cache the Silverman bandwidth for a reference sample
#'
#' Wraps \code{\link{kde_bandwidth_silverman}} with a simple in-memory cache
#' keyed on a checksum of \code{ref_sample}.  Subsequent calls with the same
#' reference matrix return the cached bandwidth matrix instantly, which is
#' particularly useful during likelihood optimisation where the background
#' sample is fixed across hundreds of function evaluations.
#'
#' @param ref_sample Numeric matrix or data frame of reference (background)
#'   points (rows = observations, columns = variables).
#' @return Bandwidth matrix H (p x p) suitable for passing to
#'   \code{\link{kde_eval_cached}}.
#' @export
#'
#' @examples
#' H <- kde_bandwidth_cached(samMPts)
kde_bandwidth_cached <- function(ref_sample) {
  ref_sample <- as.matrix(ref_sample)
  # Collision-resistant checksum: combine dimensions, column sums, row sums
  # of the first and last rows, and the Frobenius norm (sum of squares).
  # This avoids adding a `digest` dependency while being highly unlikely
  # to produce false cache hits for distinct ecological background samples.
  key <- paste0(
    nrow(ref_sample), "x", ncol(ref_sample), "_",
    paste(round(colSums(ref_sample), 8), collapse = "_"), "_",
    paste(round(ref_sample[1L, ], 8), collapse = "_"), "_",
    paste(round(ref_sample[nrow(ref_sample), ], 8), collapse = "_"), "_",
    round(sum(ref_sample^2), 8)
  )
  if (!exists(key, envir = .kde_cache, inherits = FALSE)) {
    .kde_cache[[key]] <- kde_bandwidth_silverman(ref_sample)
  }
  .kde_cache[[key]]
}

#' Invalidate the KDE bandwidth cache
#'
#' Removes all cached bandwidth matrices.  Call this if you switch to a
#' different background sample within the same R session.
#'
#' @return \code{NULL} invisibly.
#' @export
#'
#' @examples
#' kde_cache_clear()
kde_cache_clear <- function() {
  rm(list = ls(.kde_cache), envir = .kde_cache)
  invisible(NULL)
}

#' Evaluate multivariate Gaussian KDE with cached bandwidth
#'
#' A convenience wrapper around \code{\link{kde_eval_cpp}} that automatically
#' computes (and caches) the Silverman bandwidth for \code{ref_sample} via
#' \code{\link{kde_bandwidth_cached}}.  Pass a pre-computed \code{H} to
#' bypass the cache and use a custom bandwidth.
#'
#' The computation over query points is parallelised via OpenMP when the
#' package was compiled with OpenMP support.
#'
#' @param query Numeric matrix or data frame of evaluation points (m x p).
#' @param ref_sample Numeric matrix or data frame of reference points (n x p).
#' @param H Optional bandwidth matrix (p x p).  When \code{NULL} (default)
#'   the Silverman rule is applied and the result is cached.
#' @return Numeric vector of length m containing log-density values.
#' @export
#'
#' @examples
#' log_dens <- kde_eval_cached(spOccPnts, samMPts)
kde_eval_cached <- function(query, ref_sample, H = NULL) {
  query      <- as.matrix(query)
  ref_sample <- as.matrix(ref_sample)
  if (is.null(H)) {
    H <- kde_bandwidth_cached(ref_sample)
  }
  kde_eval_cpp(query, ref_sample, H)
}
