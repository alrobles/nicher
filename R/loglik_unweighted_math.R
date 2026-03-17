#' Log-likelihood for the unweighted Gaussian niche model (R math implementation)
#'
#' Computes the Gaussian log-likelihood for the unweighted niche model given
#' presence points, background points, and Cholesky-factored covariance parameters.
#'
#' @param sam1 Data frame or matrix of presence points (rows = observations, columns = variables).
#' @param sam2 Data frame or matrix of background environmental points.
#' @param mu Numeric vector of means (optimum) of length p.
#' @param L Lower triangular matrix from the Cholesky decomposition of the covariance matrix S
#'          (S = L \%*\% t(L)). Diagonal entries must be positive.
#' @return Log-likelihood (scalar).
#' @export
#'
#' @examples
#' par <- get_ellip_par(spOccPnts)
#' L <- t(chol(par$S))
#' loglik_unweighted_math(spOccPnts, samMPts, par$mu, L)
loglik_unweighted_math <- function(sam1, sam2, mu, L) {
  sam1 <- as.matrix(sam1)
  sam2 <- as.matrix(sam2)

  S <- L %*% t(L)
  logdet <- as.numeric(determinant(S, logarithm = TRUE)$modulus)

  q1 <- stats::mahalanobis(sam1, center = mu, cov = S, inverted = FALSE)
  q2 <- stats::mahalanobis(sam2, center = mu, cov = S, inverted = FALSE)

  n1 <- length(q1)
  n2 <- length(q2)

  logL <- -0.5 * (n1 * logdet + sum(q1)) + 0.5 * (n2 * logdet + sum(q2))

  return(logL)
}
