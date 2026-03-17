#' Log-likelihood for the unweighted Gaussian niche model
#'
#' Computes the Gaussian log-likelihood for the unweighted niche model given
#' presence points, background points, and Cholesky-factored covariance parameters.
#' Delegates to the C++ engine \code{\link{loglik_unweighted_cpp}}.
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
  loglik_unweighted_cpp(
    sam1 = as.matrix(sam1),
    sam2 = as.matrix(sam2),
    mu   = mu,
    L    = L
  )
}
