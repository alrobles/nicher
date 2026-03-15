#' Negative Log-Likelihood for the Weighted-Normal Niche Model (Math Scale)
#'
#' Computes the negative log-likelihood of the *weighted-normal* ecological
#' niche model described by Jiménez & Soberón (2022), using the parameterization
#' on the unconstrained “math” scale. The parameter vector \code{theta} contains:
#' \itemize{
#'   \item First \code{p} elements: \code{mu} (centroid)
#'   \item Next \code{p} elements: \code{log_sigma} (log of standard deviations)
#'   \item Last \code{p(p-1)/2} elements: parameters \code{v} for the correlation
#'         matrix via a C-vine LKJ transformation
#' }
#'
#' The covariance matrix \code{Sigma} is constructed as
#' \deqn{
#'   \Sigma = \mathrm{diag}(\sigma)\, R\, \mathrm{diag}(\sigma),
#' }
#' where \eqn{\sigma = \exp(\log \sigma)}, and \eqn{R = LL^\top} with
#' \eqn{L = \text{cvine\_cholesky}(v)}.
#'
#' The weighted-normal log-likelihood is:
#' \deqn{
#'   -\log L
#'   = \frac{1}{2}\sum_i q_1(x_i)
#'     - \sum_i \log w(x_i)
#'     + n\,\log\!\left(\sum_j w(y_j)\, e^{-q_2(y_j)/2}\right),
#' }
#' where the weights \eqn{w(\cdot)} are KDEs computed from \code{env_m}. The
#' Mahalanobis terms \eqn{q_1} and \eqn{q_2} correspond to presences and M-sample
#' points, respectively.
#'
#' @param theta Numeric vector of length \code{2p + p(p-1)/2} containing
#'   \code{mu}, \code{log_sigma}, and \code{v} in sequence.
#' @param env_occ Data frame with environmental values at presence points
#'   (size \code{n_occ x p}).
#' @param env_m Data frame with environmental values from the M area
#'   (size \code{n_m x p}); used both for KDE and the denominator.
#' @param eta Numeric. LKJ shape parameter for the C-vine transformation
#'   (default \code{1}, gives uniform over correlations).
#' @param neg Logical. If \code{TRUE} (default), returns the negative
#'   log-likelihood (for minimization).
#'
#' @return Negative log-likelihood value (scalar).
#' @export
#'
#' @examples
#' theta <- start_theta(example_env_occ_2d)
#' loglik_niche_math_weighted(theta, example_env_occ_2d, example_env_m_2d)
loglik_niche_math_weighted <- function(theta, env_occ, env_m, eta = 1, neg = TRUE) {
  
  # Dimension
  p <- ncol(env_occ)
  
  # Extract components
  mu        <- theta[1:p]
  log_sigma <- theta[(p + 1):(2 * p)]
  sigma     <- exp(log_sigma)
  v         <- if (p > 1) theta[(2 * p + 1):length(theta)] else numeric(0)
  
  # Build correlation matrix via C-vine Cholesky factor
  # (cvine_cholesky returns lower-triangular L such that R = L %*% t(L))
  l_mat <- cvine_cholesky(v, d = p, eta = eta)
  r_mat <- tcrossprod(l_mat) # correlation matrix
  
  # Construct covariance matrix
  # Equivalent to diag(sigma) %*% R %*% diag(sigma) but faster
  s_mat <- r_mat * outer(sigma, sigma)
  
  # Call weighted log-likelihood on the natural scale
  loglik_niche_weighted(
    mu     = mu,
    s_mat  = s_mat,
    env_occ = env_occ,
    env_m   = env_m,
    neg     = neg
  )
}
