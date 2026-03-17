#' Log-likelihood for the unweighted Gaussian niche model (biological-scale)
#'
#' User-facing wrapper that accepts biological-scale parameters directly
#' (niche centroid \code{mu} and lower-triangular Cholesky factor \code{L}) and
#' delegates to \code{\link{loglik_unweighted_math}}.  This is the complement of
#' \code{\link{loglik_math_niche}}, which accepts parameters in the flattened
#' math-scale (log-Cholesky) parameterisation used by optimisers.
#'
#' @param sam1 Data frame or matrix of presence points
#'   (rows = observations, columns = variables).
#' @param sam2 Data frame or matrix of background environmental points.
#' @param mu Numeric vector of means (niche optimum) of length \code{p}.
#' @param L Lower-triangular Cholesky factor of the covariance matrix
#'   (\eqn{\Sigma = L L^\top}{S = L \%*\% t(L)}).  Diagonal entries must be
#'   positive; this is validated implicitly by the downstream
#'   \code{\link{loglik_unweighted_math}} call.
#' @return Log-likelihood (scalar double).
#' @export
#'
#' @examples
#' par <- get_ellip_par(spOccPnts)
#' L <- t(chol(par$S))
#' loglik_bio_niche(spOccPnts, samMPts, par$mu, L)
loglik_bio_niche <- function(sam1, sam2, mu, L) {
  loglik_unweighted_math(sam1, sam2, mu, L)
}
