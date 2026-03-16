#' Negative log-likelihood (presence-only, math scale)
#'
#' Computes the negative log-likelihood of a multivariate normal model for presence-only data.
#' This is the simplest form, without background correction (Jiménez et al. 2019 style but without M).
#' Useful for quick estimates or comparisons.
#'
#' @param theta Numeric vector of unconstrained parameters (math scale).
#' @param env_occ Data frame with environmental values at presence points.
#' @param eta Numeric, shape parameter for LKJ‑C‑vine prior (default 1).
#' @param neg Logical, return negative log-likelihood? (default TRUE)
#'
#' @return Scalar numeric: (negative) log-likelihood.
#' @export
#'
#' @examples
#' \dontrun{
#' theta <- start_theta(example_env_occ_2d)
#' ll <- loglik_niche_math_presence_only(theta, example_env_occ_2d)
#' }
loglik_niche_math_presence_only <- function(theta, env_occ, eta = 1, neg = TRUE) {
  p <- ncol(env_occ)
  mu <- theta[1:p]
  log_sigma <- theta[(p + 1):(2 * p)]
  sigma <- exp(log_sigma)
  v <- if (p > 1) theta[(2 * p + 1):length(theta)] else numeric(0)
  
  L_corr <- cvine_cholesky(v, d = p, eta = eta)
  L_cov <- diag(sigma) %*% L_corr
  
  env_occ <- as.matrix(env_occ)
  
  val <- loglik_niche_presence_only_cpp(mu, L_cov, env_occ)
  if (neg) val else -val
}