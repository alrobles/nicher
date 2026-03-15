#' Negative log-likelihood (weighted-normal, Jiménez & Soberón 2022)
#'
#' Implements the weighted likelihood:
#'   -log L = 0.5*sum(q1) - sum(log w_occ) + n * log( sum_j w_m[j] * exp(-0.5*q2_j) )
#' where q1 uses env_occ, q2 uses env_m, and w(.) is a KDE over env_m.
#'
#' @param mu   Numeric length-p vector (centroid).
#' @param s_mat Numeric p x p covariance matrix (positive-definite).
#' @param env_occ Data frame/matrix (n_occ x p): environmental values at presences.
#' @param env_m   Data frame/matrix (n_m   x p): environmental values from M.
#' @param neg    Logical, return negative log-likelihood (default TRUE).
#' @return Numeric scalar (neg log-likelihood if neg=TRUE, else log-likelihood).
#' @export
#' @examples
#' 
#' pars <- get_ellipsoid_pars(example_env_occ_2d)  # lista(mu, s_mat)
#' negll_w <- loglik_niche_weighted(
#'   mu = pars$mu,
#'   s_mat = pars$s_mat,
#'   env_occ = example_env_occ_2d,
#'   env_m   = example_env_m_2d
#'   )
loglik_niche_weighted <- function(mu, s_mat, env_occ, env_m, neg = TRUE) {
  env_occ <- as.matrix(env_occ)
  env_m   <- as.matrix(env_m)
  
  # Distancias de Mahalanobis
  q1 <- stats::mahalanobis(env_occ, center = mu, cov = s_mat, inverted = FALSE)
  q2 <- stats::mahalanobis(env_m,   center = mu, cov = s_mat, inverted = FALSE)
  
  # KDE en M (siempre base R)
  w_m   <- kde_gaussian(env_m,   env_m)
  w_occ <- kde_gaussian(env_occ, env_m)
  
  # log-sum-exp para estabilidad
  log_sum_exp <- function(a) {
    m <- max(a)
    m + log(sum(exp(a - m)))
  }
  
  n <- length(q1)
  num <- 0.5 * sum(q1) - sum(log(w_occ))
  den <- n * log_sum_exp(log(w_m) - 0.5 * q2)
  neg_log <- num + den
  
  if (neg) neg_log else -neg_log
}
