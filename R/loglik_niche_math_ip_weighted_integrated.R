#' Negative log-likelihood (inverse-probability-weighted normal, integrated C++)
#'
#' Low-level C++ bridge for the inverse-probability-weighted (IPW) normal
#' niche model. Computes KDE weights and Mahalanobis distances entirely
#' in C++ to minimise R overhead. This is the workhorse called by
#' \code{\link{loglik_niche_math_ip_weighted}}.
#'
#' The objective function is described in detail in
#' \code{\link{loglik_niche_math_ip_weighted}}.
#'
#' This function was formerly called
#' \code{loglik_niche_math_kde_bias_corrected_integrated}.
#'
#' @param theta Numeric vector of unconstrained parameters (math scale):
#'   \eqn{[\mu, \log\sigma, v]} of length \eqn{2p + p(p-1)/2}.
#' @param env_occ Data frame or matrix (\eqn{n \times p}) of environmental
#'   values at presence points.
#' @param env_m Data frame or matrix (\eqn{n_m \times p}) of background
#'   environmental values from the accessible area M.
#' @param eta Numeric scalar, shape parameter for the LKJ C-vine prior on the
#'   correlation matrix (default 1 = uniform over correlations).
#' @param neg Logical. If \code{TRUE} (default), returns the negative
#'   log-likelihood (suitable for minimisation).
#' @param den_idx Integer vector of 1-based row indices for the denominator
#'   subsample. If \code{NULL} (default), uses all rows of \code{env_m}.
#' @param kde_idx Integer vector of 1-based row indices for the KDE reference
#'   subsample. If \code{NULL}, uses all rows of \code{env_m}.
#' @param precomp_w_den Numeric vector of precomputed KDE weights for the
#'   denominator points (must match \code{den_idx} in length). If provided,
#'   KDE for the denominator is not recomputed.
#'
#' @return Scalar numeric: the (negative) log-likelihood.
#'
#' @seealso [loglik_niche_math_ip_weighted()] for the user-facing
#'   wrapper with automatic subsampling.
#'
#' @examples
#' \donttest{
#' theta <- start_theta(example_env_occ_2d)
#' ll <- loglik_niche_math_ip_weighted_integrated(
#'   theta   = theta,
#'   env_occ = example_env_occ_2d,
#'   env_m   = example_env_m_2d,
#'   eta     = 1,
#'   neg     = TRUE
#' )
#' print(ll)
#' }
#'
#' @export
loglik_niche_math_ip_weighted_integrated <- function(theta, env_occ, env_m, eta = 1, neg = TRUE,
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

  # Convert to matrix
  env_occ <- as.matrix(env_occ)
  env_m <- as.matrix(env_m)

  # Coerce indices to integer to prevent Rcpp type mismatch errors
  if (!is.null(den_idx)) den_idx <- as.integer(den_idx)
  if (!is.null(kde_idx)) kde_idx <- as.integer(kde_idx)
  if (!is.null(precomp_w_den)) precomp_w_den <- as.numeric(precomp_w_den)

  # Call C++ with the new argument
  loglik_niche_ip_weighted_cpp(mu, L_cov, env_occ, env_m, den_idx, kde_idx, precomp_w_den, neg)
}
