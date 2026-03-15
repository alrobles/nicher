#' Negative log-likelihood (weighted-normal, Jiménez & Soberón 2022) with MC subsampling
#'
#' Implements the weighted likelihood:
#'   -log L = 0.5*sum(q1) - sum(log w_occ) + n * log( sum_j w_m[j] * exp(-0.5*q2_j) )
#' where q1 uses env_occ, q2 uses a (possibly subsampled) set from env_m, and
#' w(.) is a KDE over env_m (optionally subsampled for KDE construction).
#'
#' @param mu   Numeric length-p vector (centroid).
#' @param s_mat Numeric p x p covariance matrix (positive-definite).
#' @param env_occ Data frame/matrix (n_occ x p): environmental values at presences.
#' @param env_m   Data frame/matrix (n_m   x p): environmental values from M.
#' @param m_subsample Optional integer or fraction. If NULL (default), uses all rows of
#'   env_m in the denominator. If 0 < value < 1, uses `floor(value * nrow(env_m))`
#'   rows (without replacement). If >= 1, uses `min(value, nrow(env_m))` rows.
#' @param m_kde_subsample Optional integer or fraction controlling the KDE reference
#'   sample from env_m. Same rules as `m_subsample`. By default (NULL), KDE is built
#'   using all rows of env_m.
#' @param seed Optional integer seed for reproducible subsampling (set only here).
#' @param neg  Logical, return negative log-likelihood (default TRUE).
#'
#' @return Numeric scalar (neg log-likelihood if neg=TRUE, else log-likelihood).
#' @examples
#' pars <- get_ellipsoid_pars(example_env_occ_2d)
#' # Full M (no subsampling)
#' loglik_niche_weighted(pars$mu, pars$s_mat, example_env_occ_2d, example_env_m_2d)
#' # Subsample 2,000 puntos de M para el denominador, y 5,000 para el KDE
#' loglik_niche_weighted(pars$mu, pars$s_mat, example_env_occ_2d, example_env_m_2d,
#'                       m_subsample = 2000, m_kde_subsample = 5000, seed = 123)
#' @export
loglik_niche_weighted <- function(mu, s_mat, env_occ, env_m,
                                  m_subsample = NULL,
                                  m_kde_subsample = NULL,
                                  seed = NULL,
                                  neg = TRUE) {
  
  # --- prepare matrices -------------------------------------------------------
  env_occ <- as.matrix(env_occ)
  env_m   <- as.matrix(env_m)
  if (ncol(env_occ) != ncol(env_m)) {
    stop("env_occ and env_m must have the same number of columns (same variables).")
  }
  
  # --- Mahalanobis terms ------------------------------------------------------
  q1 <- stats::mahalanobis(env_occ, center = mu, cov = s_mat, inverted = FALSE)
  
  # Denominator set (M_den): possibly subsampled from env_m
  n_m <- nrow(env_m)
  if (!is.null(seed)) {
    old_seed <- .Random.seed; on.exit({ if (exists("old_seed")) .Random.seed <<- old_seed }, add = TRUE)
    set.seed(seed)
  }
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
    M_den   <- env_m[idx_den, , drop = FALSE]
  } else {
    M_den   <- env_m
  }
  q2 <- stats::mahalanobis(M_den, center = mu, cov = s_mat, inverted = FALSE)
  
  # KDE reference set (M_kde): possibly subsampled from env_m (default: all env_m)
  n_kde <- pick_size(m_kde_subsample, n_m)
  if (n_kde < n_m) {
    idx_kde <- sample.int(n_m, size = n_kde, replace = FALSE)
    M_kde   <- env_m[idx_kde, , drop = FALSE]
  } else {
    M_kde   <- env_m
  }
  
  # --- Weights via KDE (C++ Scott) -------------------------------------------
  # w_m: evaluate KDE at denominator points; w_occ: evaluate at presences
  w_m   <- kde_gaussian(M_den,   M_kde)
  w_occ <- kde_gaussian(env_occ, M_kde)
  
  # --- log-sum-exp helper -----------------------------------------------------
  log_sum_exp <- function(a) {
    m <- max(a); m + log(sum(exp(a - m)))
  }
  
  # --- Weighted likelihood ----------------------------------------------------
  n   <- length(q1)
  num <- 0.5 * sum(q1) - sum(log(w_occ))
  den <- n * log_sum_exp(log(w_m) - 0.5 * q2)
  neg_log <- num + den
  
  if (neg) neg_log else -neg_log
}
