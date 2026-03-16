#' Negative log-likelihood (weighted-normal, math scale)
#'
#' Computes the weighted-normal log-likelihood (Jiménez & Soberón 2022) using
#' an efficient C++ implementation with integrated KDE. This function is a
#' wrapper that converts subsampling parameters to indices and calls the
#' internal integrated version.
#'
#' @param theta Numeric vector of unconstrained parameters (math scale).
#' @param env_occ Data frame with environmental values at presence points.
#' @param env_m Data frame with environmental values from the M area.
#' @param eta Numeric, shape parameter for LKJ‑C‑vine prior (default 1).
#' @param neg Logical, return negative log-likelihood? (default TRUE)
#' @param m_subsample Optional integer or fraction for denominator subsample.
#' @param m_kde_subsample Optional integer or fraction for KDE reference subsample.
#' @param seed Optional integer seed for reproducible subsampling.
#' @param ... Additional arguments (not used, for compatibility).
#'
#' @return Scalar numeric: (negative) log-likelihood.
#' @export
#'
#' @examples
#' \dontrun{
#' theta <- start_theta(example_env_occ_2d)
#' loglik_niche_math_weighted(theta, example_env_occ_2d, example_env_m_2d,
#'                            m_subsample = 2000, m_kde_subsample = 5000, seed = 123)
#' }
loglik_niche_math_weighted <- function(theta, env_occ, env_m, eta = 1, neg = TRUE,
                                       m_subsample = NULL, m_kde_subsample = NULL,
                                       seed = NULL, ...) {
  # Dimensiones y validación básica
  p <- ncol(env_occ)
  if (p != ncol(env_m)) stop("env_occ and env_m must have the same number of columns")
  
  # Generar índices de submuestreo si se solicitan
  n_m <- nrow(env_m)
  if (!is.null(seed)) set.seed(seed)
  
  pick_size <- function(x, nmax) {
    if (is.null(x)) return(NULL)
    if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x <= 0) {
      stop("m_subsample/m_kde_subsample must be a single positive numeric (fraction or count).")
    }
    if (x < 1) max(1L, floor(x * nmax)) else min(nmax, as.integer(round(x)))
  }
  
  n_den <- pick_size(m_subsample, n_m)
  den_idx <- if (!is.null(n_den) && n_den < n_m) sample.int(n_m, n_den) else NULL
  
  n_kde <- pick_size(m_kde_subsample, n_m)
  kde_idx <- if (!is.null(n_kde) && n_kde < n_m) sample.int(n_m, n_kde) else NULL
  
  # Llamar a la versión integrada
  loglik_niche_math_weighted_integrated(theta, env_occ, env_m, eta, neg,
                                        den_idx, kde_idx)
}