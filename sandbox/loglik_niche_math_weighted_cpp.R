#' Negative log-likelihood (weighted-normal, math scale, Cholesky version)
#'
#' Computes the weighted-normal log-likelihood (Jiménez & Soberón 2022) using
#' the unconstrained ("math") parameterization and passing the Cholesky factor
#' directly to C++ for efficient calculation. Handles optional subsampling of
#' the M background for the denominator and for KDE estimation.
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
#' @param eta Numeric, shape parameter for the LKJ‑C‑vine prior (default 1).
#' @param neg Logical. If `TRUE` (default) returns the negative log‑likelihood
#'            (suitable for minimization); if `FALSE` returns the positive log‑likelihood.
#' @param m_subsample Optional integer or fraction for subsampling the denominator
#'   from `env_m`. If `NULL` (default), uses all rows. If `0 < value < 1`, uses
#'   `floor(value * nrow(env_m))` rows; if `>= 1`, uses `min(value, nrow(env_m))` rows.
#' @param m_kde_subsample Optional integer or fraction for subsampling the KDE reference
#'   set from `env_m`. Same rules as `m_subsample`.
#' @param seed Optional integer seed for reproducible subsampling.
#'
#' @return Scalar numeric value: the (negative) log-likelihood.
#'
#' @examples
#' \dontrun{
#' theta <- start_theta(example_env_occ_2d)
#' # Weighted log-likelihood (math scale) using Cholesky version
#' ll <- loglik_niche_math_weighted_cpp(theta,
#'                                      env_occ = example_env_occ_2d,
#'                                      env_m   = example_env_m_2d,
#'                                      m_subsample = 2000,
#'                                      m_kde_subsample = 5000,
#'                                      seed = 123)
#' print(ll)
#' }
#' @export
loglik_niche_math_weighted_cpp <- function(theta, env_occ, env_m, eta = 1, neg = TRUE,
                                           m_subsample = NULL, m_kde_subsample = NULL,
                                           seed = NULL) {
  # Dimensiones y validación básica
  p <- ncol(env_occ)
  if (p != ncol(env_m)) stop("env_occ and env_m must have the same number of columns")
  if (length(theta) != 2*p + p*(p-1)/2) {
    stop("theta has incorrect length for p = ", p)
  }
  
  # Extraer componentes
  mu <- theta[1:p]
  log_sigma <- theta[(p + 1):(2 * p)]
  sigma <- exp(log_sigma)
  v <- if (p > 1) theta[(2 * p + 1):length(theta)] else numeric(0)
  
  # Construir factor de Cholesky de la correlación y luego de covarianza
  L_corr <- cvine_cholesky(v, d = p, eta = eta)
  L_cov <- diag(sigma) %*% L_corr
  
  # Subsampling de M (código adaptado de loglik_niche_weighted)
  n_m <- nrow(env_m)
  if (!is.null(seed)) set.seed(seed)
  
  pick_size <- function(x, nmax) {
    if (is.null(x)) return(nmax)
    if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x <= 0) {
      stop("m_subsample/m_kde_subsample must be a single positive numeric (fraction or count).")
    }
    if (x < 1) max(1L, floor(x * nmax)) else min(nmax, as.integer(round(x)))
  }
  
  n_den <- pick_size(m_subsample, n_m)
  if (n_den < n_m) {
    idx_den <- sample.int(n_m, size = n_den, replace = FALSE)
    M_den <- env_m[idx_den, , drop = FALSE]
  } else {
    M_den <- env_m
  }
  
  n_kde <- pick_size(m_kde_subsample, n_m)
  if (n_kde < n_m) {
    idx_kde <- sample.int(n_m, size = n_kde, replace = FALSE)
    M_kde <- env_m[idx_kde, , drop = FALSE]
  } else {
    M_kde <- env_m
  }
  
  # Calcular pesos KDE (en R, usando la función existente)
  w_occ <- kde_gaussian(env_occ, M_kde)
  w_den <- kde_gaussian(M_den, M_kde)
  
  # Convertir a matrices (para evitar problemas de tipos en C++)
  env_occ_mat <- as.matrix(env_occ)
  M_den_mat <- as.matrix(M_den)
  
  # Llamar a la función C++
  val <- loglik_niche_weighted_chol_cpp(mu, L_cov, env_occ_mat, M_den_mat, w_occ, w_den)
  
  if (neg) val else -val
}