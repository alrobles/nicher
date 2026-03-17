#' Negative log-likelihood for niche model on unconstrained (math) scale
#'
#' @param theta Numeric vector of unconstrained parameters:
#'   \itemize{
#'     \item First `p` elements: `mu` (centroid)
#'     \item Next `p` elements: `log_sigma` (log of standard deviations)
#'     \item Last `p*(p-1)/2` elements: parameters for correlation matrix via C-vine
#'   }
#' @param env_occ Data frame with environmental values at presence points.
#' @param env_m Data frame with environmental values from the M area.
#' @param eta Numeric. LKJ shape parameter for the C-vine transformation
#' (default 1, gives uniform over correlations).
#' @param neg Logical. If `TRUE` (default), returns negative log-likelihood
#' (for minimization).
#'
#' @return Negative log-likelihood value.
#' @export
#' @examples
#' theta <- start_theta(example_env_occ_2d)
#' loglik_niche_math(theta, example_env_occ_2d, example_env_m_2d)
loglik_niche_math <- function(theta, env_occ, env_m, eta = 1, neg = TRUE) {
  
  p <- ncol(env_occ)
  
  # Extract components
  mu <- theta[1:p]
  log_sigma <- theta[(p + 1):(2 * p)]
  sigma <- exp(log_sigma)
  v <- theta[(2 * p + 1):length(theta)]
  
  # Build correlation matrix via C-vine Cholesky factor
  # (cvine_cholesky returns lower-triangular L such that R = L %*% t(L))
  l_mat <- cvine_cholesky(v, d = p, eta = eta)
  r_mat <- tcrossprod(l_mat)   # correlation matrix
  
  # Construct covariance matrix
  # Equivalent to diag(sigma) %*% R %*% diag(sigma) but faster
  s_mat <- r_mat * outer(sigma, sigma)
  
  # Call original log-likelihood
  loglik_niche(mu = mu, s_mat = s_mat, env_occ = env_occ, env_m = env_m, neg = neg)
}