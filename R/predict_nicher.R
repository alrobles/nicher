# S3 predict method: fitted nicher -> habitat-suitability raster.

#' Habitat-suitability raster from a fitted \code{nicher} object
#'
#' \code{predict} method for objects of class \code{"nicher"} returned
#' by \code{\link{optimize_niche}}. Reconstructs the niche centroid
#' \eqn{\mu} and covariance \eqn{\Sigma} from the optimizer's
#' \code{best$theta} (mu, log_sigma, v) parameterization, then evaluates
#' the standardized Gaussian suitability map of Jimenez et al. (2022,
#' Eq. 2) over an environmental \code{\link[terra]{SpatRaster}}.
#'
#' The same Gaussian suitability surface is used for every fitted family:
#' \code{"weighted"}, \code{"ip_weighted"}, \code{"presence_only"},
#' and the skew-normal / skew-t variants all project the fitted
#' \eqn{(\mu, \Sigma)} geometry. For skew fits, the skewness and tail
#' parameters affect fitting but are not used by this Gaussian projection.
#'
#' @param object A \code{"nicher"} object returned by
#'   \code{\link{optimize_niche}}.
#' @param env A multi-layer \code{\link[terra]{SpatRaster}}, one layer
#'   per environmental variable. When \code{object$var_names} is
#'   non-\code{NULL} (the default for fits produced by
#'   \code{optimize_niche()} on a column-named \code{env_occ}), the
#'   layers of \code{env} are reordered to match \code{var_names} and
#'   any name not present in \code{env} produces an error. When
#'   \code{object$var_names} is \code{NULL} the layers are used in
#'   their existing order. This lets you call \code{predict} on a
#'   future-climate \code{SpatRaster} that contains the same variables
#'   in any order, including extra layers.
#' @param ... Currently ignored.
#' @param return_log Logical. \code{FALSE} (default) returns suitability
#'   in \eqn{(0, 1]}; \code{TRUE} returns log-suitability in
#'   \eqn{(-\infty, 0]}.
#' @param threads Integer. Threads for the inner C++ kernel. Default
#'   \code{RcppParallel::defaultNumThreads()}.
#' @param output Character. File path for the output GeoTIFF. The
#'   empty string \code{""} (default) returns an in-memory raster.
#' @param overwrite Logical. Forwarded to \code{\link{habitat_suitability}}.
#' @param wopt List. Forwarded to \code{\link{habitat_suitability}}.
#'
#' @return A one-layer \code{SpatRaster} named \code{"suitability"}
#'   (or \code{"log_suitability"}). Returned invisibly when
#'   \code{output} is non-empty.
#'
#' @details
#' Internally:
#' \enumerate{
#'   \item \eqn{\mu = \theta_{1:p}}.
#'   \item \eqn{\sigma = \exp(\theta_{p+1:2p})}.
#'   \item C-vine partial-correlation block
#'     \eqn{v = \theta_{2p+1:\,length(\theta)}} is mapped to a
#'     correlation Cholesky factor via
#'     \code{\link{cvine_cholesky}}.
#'   \item \eqn{L = \mathrm{diag}(\sigma)\, L_{corr}},
#'     \eqn{\Sigma = L L^\top}.
#' }
#'
#' @seealso \code{\link{habitat_suitability}},
#'   \code{\link{niche_suitability_cpp}},
#'   \code{\link{cvine_cholesky}}.
#'
#' @examples
#' \donttest{
#'   fit <- optimize_niche(
#'     env_occ    = example_env_occ_2d,
#'     env_m      = example_env_m_2d,
#'     num_starts = 5L,
#'     likelihood = "weighted"
#'   )
#'   r1 <- terra::rast(nrows = 10, ncols = 10)
#'   r2 <- terra::rast(nrows = 10, ncols = 10)
#'   terra::values(r1) <- rnorm(100)
#'   terra::values(r2) <- rnorm(100)
#'   env <- c(r1, r2)
#'   names(env) <- colnames(example_env_occ_2d)
#'   suit <- predict(fit, env)
#' }
#'
#' @export
predict.nicher <- function(
    object,
    env,
    ...,
    return_log = FALSE,
    threads    = RcppParallel::defaultNumThreads(),
    output     = "",
    overwrite  = FALSE,
    wopt       = list()
) {
  if (!inherits(object, "nicher")) {
    stop("`object` must be a 'nicher' object returned by optimize_niche().")
  }
  if (!inherits(env, "SpatRaster")) {
    stop("`env` must be a terra SpatRaster (one layer per variable).")
  }

  theta <- object$best$theta
  if (!is.numeric(theta) || length(theta) < 2L) {
    stop("object$best$theta is missing or malformed.")
  }

  pars <- .recover_mu_sigma(object)
  p <- pars$p

  # Reorder env layers by name when names are available on the fit.
  env <- .reorder_env_for_predict(env, object$var_names, p)

  habitat_suitability(
    param      = list(mu = pars$mu, Sigma = pars$Sigma),
    env        = env,
    output     = output,
    overwrite  = overwrite,
    return_log = return_log,
    threads    = threads,
    wopt       = wopt
  )
}

# ---------------------------------------------------------------------------
# Internal: reorder a SpatRaster by variable name when names are stored on
# the fit. With no var_names, fall back to positional ordering and check
# layer count.
# ---------------------------------------------------------------------------
.reorder_env_for_predict <- function(env, var_names, p) {
  if (is.null(var_names)) {
    if (terra::nlyr(env) != p) {
      stop("env has ", terra::nlyr(env),
           " layers but the fitted model expects ", p, ".")
    }
    return(env)
  }

  if (length(var_names) != p) {
    # Defensive: a malformed object stored a vector of the wrong length.
    stop("object$var_names has length ", length(var_names),
         " but the fitted model expects ", p, " variables.")
  }

  layer_names <- names(env)
  if (any(!nzchar(layer_names))) {
    stop("env has unnamed layers; predict.nicher requires named layers ",
         "matching object$var_names = ",
         paste(var_names, collapse = ", "), ".")
  }
  missing_names <- setdiff(var_names, layer_names)
  if (length(missing_names) > 0L) {
    stop("env is missing required layer(s): ",
         paste(missing_names, collapse = ", "), ".")
  }
  env[[var_names]]
}
