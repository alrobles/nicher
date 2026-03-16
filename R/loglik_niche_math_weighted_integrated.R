#' Negative log-likelihood (weighted-normal) integrated C++ version
#'
#' Computes the weighted-normal log-likelihood (Jiménez & Soberón 2022) entirely
#' in C++, including KDE calculation and Mahalanobis distances, to minimize
#' R overhead. Handles optional subsampling via integer indices.
#'
#' @param theta Numeric vector of unconstrained parameters (math scale).
#' @param env_occ Data frame with environmental values at presence points.
#' @param env_m   Data frame with environmental values from the M area.
#' @param eta Numeric, shape parameter for LKJ‑C‑vine prior (default 1).
#' @param neg Logical, return negative log-likelihood? (default TRUE)
#' @param den_idx Integer vector of indices (1-based) for denominator subsample.
#'   If NULL (default), uses all rows of env_m.
#' @param kde_idx Integer vector of indices for KDE reference subsample.
#'   If NULL, uses all rows of env_m.
#'
#' @return Scalar numeric: (negative) log-likelihood.
#'
#' @examples
#' \dontrun{
#' theta <- start_theta(example_env_occ_2d)
#' # Use all M
#' ll <- loglik_niche_math_weighted_integrated(theta, example_env_occ_2d, example_env_m_2d)
#' # Subsample
#' idx <- sample(1:nrow(example_env_m_2d), 2000)
#' ll2 <- loglik_niche_math_weighted_integrated(theta, example_env_occ_2d, example_env_m_2d,
#'                                              den_idx = idx, kde_idx = idx[1:500])
#' }
#' @export
loglik_niche_math_weighted_integrated <- function(theta, env_occ, env_m, eta = 1, neg = TRUE,
                                                  den_idx = NULL, kde_idx = NULL) {
  p <- ncol(env_occ)
  if (p != ncol(env_m)) stop("env_occ and env_m must have same columns")
  
  mu <- theta[1:p]
  log_sigma <- theta[(p + 1):(2 * p)]
  sigma <- exp(log_sigma)
  v <- if (p > 1) theta[(2 * p + 1):length(theta)] else numeric(0)
  
  L_corr <- cvine_cholesky(v, d = p, eta = eta)
  L_cov <- diag(sigma) %*% L_corr
  
  # Convertir a matriz
  env_occ <- as.matrix(env_occ)
  env_m   <- as.matrix(env_m)
  
  loglik_niche_weighted_integrated_cpp(mu, L_cov, env_occ, env_m, den_idx, kde_idx, neg)
}