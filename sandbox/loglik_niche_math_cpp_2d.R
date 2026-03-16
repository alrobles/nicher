#' Negative log-likelihood (niche model, math scale, 2D optimized version)
#'
#' Computes the negative log-likelihood of the ellipsoid niche model (Jiménez & Soberón 2019)
#' for the special case of **2 environmental dimensions**. This version uses a manually
#' optimized C++ implementation that avoids Eigen overhead and is intended for
#' benchmarking and performance comparison. For general dimension, use
#' [loglik_niche_math()] or [loglik_niche_math_cpp()].
#'
#' @param theta Numeric vector of unconstrained parameters:
#'   \itemize{
#'     \item First 2 elements: `mu` (centroid)
#'     \item Next 2 elements: `log_sigma` (log of standard deviations)
#'     \item Last 1 element: `v` (single C‑vine partial correlation parameter,
#'           because for p=2, `v` has length 1).
#'   }
#' @param env_occ Data frame with environmental values at presence points (size `n_occ x 2`).
#' @param env_m   Data frame with environmental values from the accessibility area M
#'                (size `n_m x 2`); used as the background sample.
#' @param eta Numeric, shape parameter for the LKJ‑C‑vine prior (default 1, uniform over
#'            correlation matrices). Passed to [cvine_cholesky()].
#' @param neg Logical. If `TRUE` (default) returns the negative log‑likelihood
#'            (suitable for minimization); if `FALSE` returns the positive log‑likelihood.
#'
#' @return A scalar numeric value: the (negative) log‑likelihood.
#'
#' @note This function is only valid when the environmental data have exactly two variables.
#'   It is provided for performance testing and educational purposes.
#'
#' @examples
#' \dontrun{
#' theta <- start_theta(example_env_occ_2d)
#' # Only works because example_env_occ_2d has 2 columns
#' ll <- loglik_niche_math_cpp_2d(theta,
#'                                env_occ = example_env_occ_2d,
#'                                env_m   = example_env_m_2d,
#'                                eta = 1, neg = TRUE)
#' print(ll)
#' }
#' @export
loglik_niche_math_cpp_2d <- function(theta, env_occ, env_m, eta = 1, neg = TRUE) {
  p <- ncol(env_occ)
  if (p != 2) stop("Solo para 2 dimensiones")
  mu <- theta[1:p]
  log_sigma <- theta[(p + 1):(2 * p)]
  sigma <- exp(log_sigma)
  v <- theta[(2 * p + 1):length(theta)]
  L_corr <- cvine_cholesky(v, d = p, eta = eta)
  L_cov <- diag(sigma) %*% L_corr
  val <- loglik_niche_chol_cpp_2d(mu, L_cov, env_occ, env_m)
  if (neg) val else -val
}