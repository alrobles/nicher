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
#' @param ... Additional arguments forwarded to \code{loglik_niche_weighted()},
#'   e.g., \code{m_subsample}, \code{m_kde_subsample}, \code{seed}.
#'
#' @return Negative log-likelihood value (scalar).
#' @export
#'
#' @examples
#' theta <- start_theta(example_env_occ_2d)
#' # Without subsampling
#' loglik_niche_math_weighted(theta, example_env_occ_2d, example_env_m_2d)
#' # With subsampling in M (Monte Carlo)
#' loglik_niche_math_weighted(theta, example_env_occ_2d, example_env_m_2d,
#'                            m_subsample = 2000, m_kde_subsample = 5000, seed = 2026)
loglik_niche_math_weighted <- function(theta, env_occ, env_m, eta = 1, neg = TRUE, ...) {
  
  # --- Basic validation ------------------------------------------------------
  if (!is.numeric(eta) || length(eta) != 1L || !(eta > 0)) {
    stop("eta must be a single positive number (> 0).")
  }
  if (ncol(env_occ) != ncol(env_m)) {
    stop("env_occ and env_m must have the same number of columns (same variables).")
  }
  
  # Dimension and theta length
  p <- ncol(env_occ)
  if (p < 1L) stop("p must be >= 1.")
  len_expected <- 2L * p + p * (p - 1L) / 2L
  if (length(theta) != len_expected) {
    stop(sprintf("theta must have length %d for p=%d (got %d).",
                 len_expected, p, length(theta)))
  }
  
  # --- Extract components ----------------------------------------------------
  mu        <- theta[1:p]
  log_sigma <- theta[(p + 1):(2 * p)]
  sigma     <- exp(log_sigma)
  v         <- if (p > 1L) theta[(2 * p + 1):length(theta)] else numeric(0)
  
  # --- Build correlation matrix via C-vine ----------------------------------
  # (cvine_cholesky returns lower-triangular L such that R = L %*% t(L))
  l_mat <- cvine_cholesky(v, d = p, eta = eta)
  r_mat <- tcrossprod(l_mat)  # correlation matrix
  
  # --- Construct covariance matrix ------------------------------------------
  # Equivalent to diag(sigma) %*% R %*% diag(sigma) but faster
  s_mat <- r_mat * outer(sigma, sigma)
  
  # --- Call weighted log-likelihood on the natural scale --------------------
  loglik_niche_weighted(
    mu      = mu,
    s_mat   = s_mat,
    env_occ = env_occ,
    env_m   = env_m,
    neg     = neg,
    ...
  )
}