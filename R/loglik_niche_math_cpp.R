#' Negative log-likelihood (niche model, math scale, Cholesky version)
#'
#' Computes the negative log-likelihood of the ellipsoid niche model (Jiménez & Soberón 2019)
#' using the unconstrained ("math") parameterization, but passes the Cholesky factor
#' directly to the C++ function `loglik_niche_chol_cpp` to avoid redundant matrix
#' decompositions. This function is intended for internal use and for comparing
#' performance with the original `loglik_niche_math`.
#'
#' @param theta Numeric vector of unconstrained parameters:
#'   \itemize{
#'     \item First `p` elements: `mu` (centroid)
#'     \item Next `p` elements: `log_sigma` (log of standard deviations)
#'     \item Last `p*(p-1)/2` elements: parameters for the correlation matrix
#'           (C‑vine partial correlations, passed to `cvine_cholesky`).
#'   }
#' @param env_occ Data frame with environmental values at presence points (size `n_occ x p`).
#' @param env_m   Data frame with environmental values from the accessibility area M
#'                (size `n_m x p`); used as the background sample.
#' @param eta Numeric, shape parameter for the LKJ‑C‑vine prior (default 1, uniform over
#'            correlation matrices). Passed to `cvine_cholesky`.
#' @param neg Logical. If `TRUE` (default) returns the negative log‑likelihood
#'            (suitable for minimization); if `FALSE` returns the positive log‑likelihood.
#' @param ... Additional arguments (ignored, for compatibility).
#'
#' @return A scalar numeric value: the (negative) log‑likelihood.
#'
#' @examples
#' \dontrun{
#' theta <- start_theta(example_env_occ_2d)
#' # Unweighted log-likelihood (math scale) using Cholesky version
#' ll <- loglik_niche_math_cpp(theta,
#'   env_occ = example_env_occ_2d,
#'   env_m = example_env_m_2d,
#'   eta = 1, neg = TRUE
#' )
#' print(ll)
#' }
#' @export
loglik_niche_math_cpp <- function(theta, env_occ, env_m, eta = 1, neg = TRUE, ...) {
  p <- ncol(env_occ)
  mu <- theta[1:p]
  log_sigma <- theta[(p + 1):(2 * p)]
  sigma <- exp(log_sigma)
  v <- if (p > 1) theta[(2 * p + 1):length(theta)] else numeric(0)

  # Build correlation Cholesky factor (lower triangular)
  L_corr <- cvine_cholesky(v, d = p, eta = eta)
  # Build covariance Cholesky factor (lower triangular)
  L_cov <- diag(sigma) %*% L_corr

  # Convert data frames to matrices to avoid Rcpp type conversion issues
  env_occ <- as.matrix(env_occ)
  env_m <- as.matrix(env_m)

  val <- loglik_niche_chol_cpp(mu, L_cov, env_occ, env_m)
  if (neg) val else -val
}
