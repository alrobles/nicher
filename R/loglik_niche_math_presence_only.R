#' Negative log-likelihood (presence-only, math scale)
#'
#' Computes the negative log-likelihood of a multivariate normal niche model
#' using only presence records, with no background correction
#' (Jimenez et al. 2019, Eq. 2 without the M-restricted denominator).
#'
#' The minimised objective is:
#' \deqn{
#'   -\log\mathcal{L} =
#'   \frac{n}{2}\log|\Sigma|
#'   + \frac{1}{2}\sum_{i=1}^{n}
#'     (\mathbf{x}_i - \mu)^\top \Sigma^{-1}(\mathbf{x}_i - \mu)
#' }
#' The \eqn{(2\pi)} normalisation constant is dropped (does not affect the
#' optimum).
#'
#' Internally, \eqn{\Sigma = L L^\top} is reconstructed from \code{theta}
#' via \code{\link{cvine_cholesky}}, and the Mahalanobis distances are
#' computed via a triangular solve \eqn{L^{-1}(\mathbf{x}_i - \mu)}.
#'
#' @param theta Numeric vector of unconstrained parameters (math scale):
#'   \eqn{[\mu, \log\sigma, v]} of length \eqn{2p + p(p-1)/2}.
#' @param env_occ Data frame or matrix (\eqn{n \times p}) of environmental
#'   values at presence points.
#' @param eta Numeric scalar, shape parameter for the LKJ C-vine prior on the
#'   correlation matrix (default 1 = uniform over correlations).
#' @param neg Logical. If \code{TRUE} (default), returns the negative
#'   log-likelihood (suitable for minimisation).
#' @param ... Additional arguments (ignored; for compatibility).
#'
#' @return Scalar numeric: the (negative) log-likelihood.
#'
#' @references
#' Jimenez, L., Soberon, J., Christen, J. A., & Soto, D. (2019).
#' On the problem of modeling a fundamental niche from occurrence data.
#' \emph{Ecological Modelling}, 397, 109823.
#'
#' @seealso [optimize_niche()] for multi-start fitting,
#'   [loglik_niche_math_cpp()] for the M-restricted model.
#'
#' @export
#'
#' @examples
#' \donttest{
#' theta <- start_theta(example_env_occ_2d)
#' ll <- loglik_niche_math_presence_only(theta, example_env_occ_2d)
#' }
loglik_niche_math_presence_only <- function(theta, env_occ, eta = 1, neg = TRUE, ...) {
  p <- ncol(env_occ)
  mu <- theta[1:p]
  log_sigma <- theta[(p + 1):(2 * p)]
  sigma <- exp(log_sigma)
  v <- if (p > 1) theta[(2 * p + 1):length(theta)] else numeric(0)

  L_corr <- cvine_cholesky(v, d = p, eta = eta)
  L_cov <- diag(sigma) %*% L_corr

  env_occ <- as.matrix(env_occ)

  val <- loglik_niche_presence_only_cpp(mu, L_cov, env_occ)
  if (neg) val else -val
}
