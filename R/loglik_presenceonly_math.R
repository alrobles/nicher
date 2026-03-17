#' Log-likelihood for the presence-only niche model (R math implementation)
#'
#' Computes the presence-only (semi-log-likelihood) for a Gaussian niche model
#' using Mahalanobis quadratic forms for presence and background points.
#'
#' @param sam1 Data frame of presence points.
#' @param sam2 Data frame of background (M) environmental points.
#' @param mu Numeric vector of means of length p.
#' @param S Covariance matrix (p x p).
#' @return Negative log-likelihood value (scalar).
#' @export
#'
#' @examples
#' par <- get_ellip_par(spOccPnts)
#' loglik_presenceonly_math(spOccPnts, samMPts, par$mu, par$S)
loglik_presenceonly_math <- function(sam1, sam2, mu, S) {
  q1 <- stats::mahalanobis(x = sam1, center = mu, cov = S, inverted = FALSE)
  q2 <- stats::mahalanobis(x = sam2, center = mu, cov = S, inverted = FALSE)
  n <- length(q1)
  0.5 * sum(q1) + n * log(sum(exp(-0.5 * q2)))
}
