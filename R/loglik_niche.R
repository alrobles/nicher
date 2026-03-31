#' Negative log likelihood of an ellipsoid corrected with environmental
#' combinations which come from the area of study (M)
#' @param mu A vector mu of parameters
#' @param s_mat The covariance matrix from environmental data frame
#' @param env_occ A data.frame containing the original sample of environmental
#' combinations that correspond to presences
#' @param env_m A data.frame containing a second random sample of environmental
#' @param neg Logical. Default TRUE, returns the negative of the
#' likelihood
#'
#' @return A negative log likelihood value
#' @export
#'
#' @examples
#' loglik_niche(
#'   mu = example_mu_vec,
#'   s_mat = example_s_mat,
#'   env_occ = example_env_occ_2d,
#'   env_m = example_env_m_2d
#' )
#' # Example with log of the likelihood
#' loglik_niche(
#'   mu = example_mu_vec,
#'   s_mat = example_s_mat,
#'   env_occ = example_env_occ_2d,
#'   env_m = example_env_m_2d,
#'   neg = FALSE
#' )
loglik_niche <- function(mu, s_mat, env_occ, env_m, neg = TRUE) {
  # quadratic terms of environment from occurrence points
  q1 <- stats::mahalanobis(
    x = env_occ,
    center = mu,
    cov = s_mat,
    inverted = FALSE
  )
  # quadratic terms of M points
  q2 <- stats::mahalanobis(
    x = env_m,
    center = mu,
    cov = s_mat,
    inverted = FALSE
  )

  # negative log-likelihood value (using log-sum-exp for numerical stability)
  n <- length(q1)
  a <- -0.5 * q2
  max_a <- max(a)
  log_sum_exp <- max_a + log(sum(exp(a - max_a)))
  neg_log <- 0.5 * sum(q1) + n * log_sum_exp
  if (neg) {
    neg_log
  } else {
    -neg_log
  }
}
