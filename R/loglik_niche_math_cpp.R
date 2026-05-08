#' Negative log-likelihood (M-restricted, math scale, Cholesky version)
#'
#' Computes the negative log-likelihood of the multivariate normal niche model
#' restricted to the set of existing environments \eqn{\mathbf{E}(t;G)}
#' (Jimenez et al. 2019, Eq. 3--5). The density at an occurrence point is
#' normalised by the sum of the density over all background points in M:
#'
#' \deqn{
#'   -\log\mathcal{L} =
#'   \frac{1}{2}\sum_{i=1}^{n}
#'     (\mathbf{x}_i - \mu)^\top \Sigma^{-1}(\mathbf{x}_i - \mu)
#'   + n \cdot \log\!\left(
#'     \sum_{\mathbf{y} \in M}
#'     \exp\!\left[-\tfrac{1}{2}
#'       (\mathbf{y} - \mu)^\top \Sigma^{-1}(\mathbf{y} - \mu)
#'     \right]
#'   \right)
#' }
#'
#' The \eqn{|\Sigma|^{-1/2}} factor cancels between numerator and denominator
#' (Eq. 3), so it does not appear in the objective. The logsumexp is computed
#' with the max-shift trick for numerical stability.
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
#' @param ... Additional arguments (ignored; for compatibility).
#'
#' @return A scalar numeric value: the (negative) log-likelihood.
#'
#' @references
#' Jimenez, L., Soberon, J., Christen, J. A., & Soto, D. (2019).
#' On the problem of modeling a fundamental niche from occurrence data.
#' \emph{Ecological Modelling}, 397, 109823.
#'
#' @seealso [optimize_niche()] for multi-start fitting,
#'   [loglik_niche_math_presence_only()] for the unconstrained (no M) version.
#'
#' @examples
#' \donttest{
#' theta <- start_theta(example_env_occ_2d)
#' ll <- loglik_niche_math_cpp(theta,
#'   env_occ = example_env_occ_2d,
#'   env_m = example_env_m_2d,
#'   eta = 1, neg = TRUE
#' )
#' print(ll)
#' }
#' @export
loglik_niche_math_cpp <- function(theta, env_occ, env_m, eta = 1, neg = TRUE, ...) {
  p <- ncol(env_occ)
  mu <- theta[1:p]
  log_sigma <- theta[(p + 1):(2 * p)]
  sigma <- exp(log_sigma)
  v <- if (p > 1) theta[(2 * p + 1):length(theta)] else numeric(0)

  # Build correlation Cholesky factor (lower triangular)
  L_corr <- cvine_cholesky(v, d = p, eta = eta)
  # Build covariance Cholesky factor (lower triangular)
  L_cov <- diag(sigma) %*% L_corr

  # Convert data frames to matrices to avoid Rcpp type conversion issues
  env_occ <- as.matrix(env_occ)
  env_m <- as.matrix(env_m)

  val <- loglik_niche_chol_cpp(mu, L_cov, env_occ, env_m)
  if (neg) val else -val
}
