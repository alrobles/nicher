#' Negative log-likelihood (presence-only model)
#'
#' Backward-compatible wrapper around \code{\link{loglik_presenceonly_math}}.
#'
#' @param mu Numeric vector of means.
#' @param S Covariance matrix.
#' @param sam1 Data frame of presence points.
#' @param sam2 Data frame of background (M) environmental points.
#' @return Negative log-likelihood value (scalar).
#' @export
#'
#' @examples
#' par <- get_ellip_par(spOccPnts)
#' negloglike_multivariable(par$mu, par$S, spOccPnts, samMPts)
negloglike_multivariable <- function(mu, S, sam1, sam2) {
  loglik_presenceonly_math(sam1, sam2, mu, S)
}
