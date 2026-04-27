# R/habitat_suitability.R
# Tiled, zero-copy habitat-suitability evaluator built on top of
# `niche_suitability_cpp` and terra's block-wise I/O.
#
# Implements the recipe documented at
# alrobles/xsdm-devel/docs/methodology/03-terra-blockwise-io.md, simplified
# for the nicher case (no time axis: each environmental variable is a
# single layer).
#
# Suitability formula (Jimenez et al. 2022, Eq. 2 + standardization
# paragraph):
#
#   S(x) = f(x; mu, Sigma) / f(mu; mu, Sigma)
#        = exp( -1/2 * (x - mu)^T Sigma^{-1} (x - mu) )    in (0, 1].

#' Tiled habitat-suitability map from an environmental \pkg{terra} stack
#'
#' Evaluates the standardized multivariate-normal density of Jimenez et
#' al. (2022, Eq. 2) at every cell of an environmental
#' \code{\link[terra]{SpatRaster}}, returning a one-layer raster of
#' suitability values in \eqn{(0, 1]} (or, optionally, log-suitability
#' in \eqn{(-\infty, 0]}).
#'
#' Memory is bounded by the largest block selected by terra's memory
#' manager: continental-scale rasters are processed without ever
#' materialising the full grid in R. The per-cell computation is run by
#' \code{\link{niche_suitability_cpp}}, an \pkg{RcppParallel} kernel
#' that takes raw pointers into terra's column-major buffer.
#'
#' @param param Named list with components
#'   \describe{
#'     \item{\code{mu}}{Numeric vector of length \code{p} — niche centroid.}
#'     \item{\code{Sigma}}{Symmetric positive-definite \code{p x p} matrix.}
#'   }
#'   The Cholesky factor of \code{Sigma} is computed once before the
#'   block loop. If \code{chol()} fails (rank-deficient \code{Sigma}),
#'   the function emits a \code{\link{warning}} and returns a raster of
#'   \code{NA}.
#' @param env A multi-layer \code{\link[terra]{SpatRaster}}: one layer
#'   per environmental variable, in the same order as \code{param$mu}.
#'   Layer names are preserved on the input but are not used by this
#'   function — variable matching is positional. Use
#'   \code{\link{predict.nicher}} when you need name-based matching.
#' @param output Character. File path for the output GeoTIFF. The empty
#'   string \code{""} (default) returns an in-memory \code{SpatRaster}.
#' @param overwrite Logical. If \code{TRUE} and \code{output} exists, it
#'   is overwritten. Default \code{FALSE}.
#' @param return_log Logical. If \code{FALSE} (default) the output is
#'   suitability \eqn{S(x) \in (0, 1]}. If \code{TRUE} the output is
#'   \eqn{\log S(x) \le 0}.
#' @param threads Integer. Number of parallel threads handed to the
#'   inner C++ kernel. Default
#'   \code{RcppParallel::defaultNumThreads()}.
#' @param wopt List. Additional write options forwarded to
#'   \code{\link[terra]{writeStart}}.
#'
#' @return A one-layer \code{SpatRaster} named \code{"suitability"} (or
#'   \code{"log_suitability"} when \code{return_log = TRUE}). Returned
#'   invisibly when \code{output} is non-empty.
#'
#' @details
#' The function follows the streaming I/O pattern from
#' \emph{xsdm-devel}'s recipe 03:
#' \enumerate{
#'   \item \code{terra::readStart} on \code{env}, paired with an
#'     \code{on.exit(terra::readStop)};
#'   \item \code{terra::writeStart} on the output, paired with an
#'     \code{on.exit(terra::writeStop)};
#'   \item for each block \code{i}: read each variable's tile via
#'     \code{terra::readValues(..., mat = TRUE)}, pack into a flat
#'     column-major numeric vector \code{(n_tile x p)}, mask
#'     \code{NA} cells, compact, call
#'     \code{niche_suitability_cpp}, scatter the result back, and
#'     \code{terra::writeValues}.
#' }
#'
#' @seealso \code{\link{niche_suitability_cpp}},
#'   \code{\link{predict.nicher}}.
#'
#' @examples
#' \dontrun{
#'   ## Hand-built (mu, Sigma):
#'   ex <- terra::rast(matrix(rnorm(100), 10), nlyrs = 2L)
#'   names(ex) <- c("bio1", "bio12")
#'   habitat_suitability(
#'     param = list(mu = c(0, 0), Sigma = diag(2)),
#'     env   = ex
#'   )
#' }
#' @export
habitat_suitability <- function(
    param,
    env,
    output      = "",
    overwrite   = FALSE,
    return_log  = FALSE,
    threads     = RcppParallel::defaultNumThreads(),
    wopt        = list()
) {
  # ----------------------------------------------------------------
  # Input validation
  # ----------------------------------------------------------------
  checkmate::assert_list(param, names = "named", any.missing = FALSE,
                          types = c("numeric", "matrix", "array"))
  checkmate::assert_subset(c("mu", "Sigma"), names(param))
  checkmate::assert_numeric(param$mu, finite = TRUE, any.missing = FALSE,
                             min.len = 1L)
  checkmate::assert_matrix(param$Sigma, mode = "numeric", any.missing = FALSE)
  if (!inherits(env, "SpatRaster")) {
    stop("`env` must be a terra SpatRaster (one layer per variable).")
  }
  checkmate::assert_string(output)
  checkmate::assert_flag(overwrite)
  checkmate::assert_flag(return_log)
  checkmate::assert_count(threads, positive = FALSE)
  checkmate::assert_list(wopt)

  mu <- as.numeric(param$mu)
  Sigma <- as.matrix(param$Sigma)
  p <- length(mu)
  if (!isTRUE(all.equal(dim(Sigma), c(p, p)))) {
    stop("param$Sigma must be a ", p, " x ", p, " matrix matching length(param$mu).")
  }
  if (terra::nlyr(env) != p) {
    stop("env must have one layer per environmental variable: ",
         "terra::nlyr(env) = ", terra::nlyr(env),
         " but length(param$mu) = ", p, ".")
  }

  # ----------------------------------------------------------------
  # Cholesky of Sigma (once). Degenerate Sigma -> warn + NA raster.
  # ----------------------------------------------------------------
  L_inv <- tryCatch({
    L <- t(chol(Sigma))                 # lower-triangular
    backsolve(L, diag(p), upper.tri = FALSE)
  }, error = function(e) {
    warning("habitat_suitability: chol(Sigma) failed (",
            conditionMessage(e),
            "); returning a raster of NA values.")
    NULL
  })

  out_name <- if (return_log) "log_suitability" else "suitability"
  out <- terra::rast(env[[1L]], nlyr = 1L)
  names(out) <- out_name

  # ----------------------------------------------------------------
  # Stream env -> kernel -> out, one block at a time.
  # ----------------------------------------------------------------
  terra::readStart(env)
  on.exit(tryCatch(terra::readStop(env), error = function(e) NULL),
          add = TRUE)

  b <- terra::writeStart(out, filename = output, overwrite = overwrite,
                         wopt = wopt)
  on.exit(tryCatch(terra::writeStop(out), error = function(e) NULL),
          add = TRUE)

  nc <- terra::ncol(env)

  for (i in seq_len(b$n)) {
    n_tile <- b$nrows[i] * nc

    if (is.null(L_inv)) {
      # Degenerate Sigma path: emit NA for the whole block.
      terra::writeValues(out, rep(NA_real_, n_tile), b$row[i], b$nrows[i])
      next
    }

    # readValues with mat=TRUE returns a (n_tile x p) matrix in
    # column-major layout — exactly what niche_suitability_cpp expects.
    tile <- terra::readValues(
      env,
      row   = b$row[i],
      nrows = b$nrows[i],
      col   = 1L,
      ncols = nc,
      mat   = TRUE
    )
    # Be defensive: terra returns whichever matrix orientation matches
    # nlyr; we want rows = pixels, cols = variables.
    if (!identical(dim(tile), c(as.integer(n_tile), as.integer(p)))) {
      tile <- matrix(as.numeric(tile), nrow = n_tile, ncol = p)
    }

    valid <- !apply(tile, 1L, anyNA)
    n_valid <- sum(valid)
    block_result <- rep(NA_real_, n_tile)

    if (n_valid > 0L) {
      if (n_valid < n_tile) {
        compact <- as.numeric(tile[valid, , drop = FALSE])
        result_valid <- niche_suitability_cpp(
          env_dat_vec  = compact,
          env_dat_dims = as.integer(c(n_valid, p)),
          mu           = mu,
          L_inv        = L_inv,
          return_log   = return_log,
          num_threads  = as.integer(threads)
        )
        block_result[valid] <- result_valid
      } else {
        block_result <- niche_suitability_cpp(
          env_dat_vec  = as.numeric(tile),
          env_dat_dims = as.integer(c(n_tile, p)),
          mu           = mu,
          L_inv        = L_inv,
          return_log   = return_log,
          num_threads  = as.integer(threads)
        )
      }
    }

    terra::writeValues(out, block_result, b$row[i], b$nrows[i])
  }

  if (nzchar(output)) invisible(out) else out
}
