# ARCHIVED LEGACY FUNCTION — not loaded by the package.
# See inst/legacy/README.md for migration guidance.
#
#' Get the negative log-likelihood value from quadratic forms of environmental
#' information from presence points and a sample from an M hypothesis
#' @param q1  quadratic terms of presence points
#' @param q2 quadratic terms of M points
#'
#' @return A single value of the negative log-likelihood
#'
#' @examples
#' par <- get_optim_par(spOccPnts)
#' q1 <- mahalanobis(spOccPnts, par$mu, par$A, inverted = TRUE)
#' q2 <- mahalanobis(samMPts, par$mu, par$A, inverted = TRUE)
#' get_negative_log(q1, q2)
get_negative_log <- function(q1, q2){
  n <- length(q1)
  0.5 * sum(q1) + n * log(sum(exp(-0.5 * q2)))
}
