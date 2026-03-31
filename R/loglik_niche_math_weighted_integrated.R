# R/loglik_niche_math_weighted_integrated.R
#' Negative log-likelihood (weighted-normal) integrated C++ version
#'
#' Computes the weighted-normal log-likelihood (Jiménez & Soberón 2022) entirely
#' in C++, including KDE calculation and Mahalanobis distances, to minimize
#' R overhead. Handles optional subsampling via integer indices. Optionally
#' accepts precomputed denominator weights to avoid recomputing KDE.
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
#' @param precomp_w_den Numeric vector of precomputed KDE weights for the
#'   denominator points (corresponding to den_idx). If provided, the function
#'   will not recompute them.
#'
#' @return Scalar numeric: (negative) log-likelihood.
#'
#' @export
loglik_niche_math_weighted_integrated <- function(theta, env_occ, env_m, eta = 1, neg = TRUE,
                                                  den_idx = NULL, kde_idx = NULL,
                                                  precomp_w_den = NULL) {
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
  env_m <- as.matrix(env_m)

  # Coerce indices to integer to prevent Rcpp type mismatch errors
  if (!is.null(den_idx)) den_idx <- as.integer(den_idx)
  if (!is.null(kde_idx)) kde_idx <- as.integer(kde_idx)
  if (!is.null(precomp_w_den)) precomp_w_den <- as.numeric(precomp_w_den)

  # Llamada a C++ con el nuevo argumento
  loglik_niche_weighted_integrated_cpp(mu, L_cov, env_occ, env_m, den_idx, kde_idx, precomp_w_den, neg)
}
