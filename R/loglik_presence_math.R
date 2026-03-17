#' Log-likelihood for the presence-only niche model
#'
#' Computes the presence-only (semi-log-likelihood) for a Gaussian niche model.
#' Delegates to the C++ engine \code{\link{loglik_presence_cpp}}.
#'
#' @param sam1 Data frame or matrix of presence points (rows = observations, columns = variables).
#' @param sam2 Data frame or matrix of background environmental points.
#' @param mu Numeric vector of means of length p.
#' @param S Covariance matrix (p x p).
#' @return Negative log-likelihood value (scalar).
#' @export
#'
#' @examples
#' par <- get_ellip_par(spOccPnts)
#' loglik_presence_math(spOccPnts, samMPts, par$mu, par$S)
loglik_presence_math <- function(sam1, sam2, mu, S) {
  loglik_presence_cpp(
    sam1 = as.matrix(sam1),
    sam2 = as.matrix(sam2),
    mu   = mu,
    S    = S
  )
}
