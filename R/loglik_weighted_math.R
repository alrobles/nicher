#' Log-likelihood for the weighted Gaussian niche model (R math implementation)
#'
#' Computes the weighted presence-only (semi-log-likelihood) for a Gaussian niche
#' model using Mahalanobis quadratic forms.  Background points can be given
#' non-uniform weights (e.g. from a Kernel Density Estimation (KDE) estimate) so that denser regions of
#' environmental space contribute proportionally more to the normalizing
#' integral.
#'
#' When all weights are equal the result reduces to
#' \code{\link{loglik_presenceonly_math}}.
#'
#' @param sam1 Data frame or matrix of presence points
#'   (rows = observations, columns = variables).
#' @param sam2 Data frame or matrix of background environmental points.
#' @param mu Numeric vector of means (optimum) of length p.
#' @param S Covariance matrix (p x p).
#' @param weights Optional numeric vector of non-negative weights for the
#'   background points in \code{sam2}.  Weights are normalized internally so
#'   that they sum to one.  When \code{NULL} (default) uniform weights are used.
#' @return Negative log-likelihood value (scalar double).
#' @export
#'
#' @examples
#' par <- get_ellip_par(spOccPnts)
#' ## Uniform weights (equivalent to loglik_presenceonly_math)
#' loglik_weighted_math(spOccPnts, samMPts, par$mu, par$S)
#'
#' ## Custom weights proportional to row index (toy example)
#' w <- seq_len(nrow(samMPts))
#' loglik_weighted_math(spOccPnts, samMPts, par$mu, par$S, weights = w)
loglik_weighted_math <- function(sam1, sam2, mu, S, weights = NULL) {
  sam2 <- as.matrix(sam2)
  n2   <- nrow(sam2)

  if (is.null(weights)) {
    weights <- rep(1.0 / n2, n2)
  } else {
    weights <- as.numeric(weights)
    if (length(weights) != n2)
      stop("'weights' must have the same length as the number of rows in sam2")
    if (any(weights < 0))
      stop("'weights' must be non-negative")
    w_sum <- sum(weights)
    if (w_sum == 0)
      stop("'weights' must not all be zero")
    weights <- weights / w_sum
  }

  q1 <- stats::mahalanobis(x = sam1, center = mu, cov = S, inverted = FALSE)
  q2 <- stats::mahalanobis(x = sam2, center = mu, cov = S, inverted = FALSE)
  n1 <- length(q1)

  # Weighted log-sum-exp for numerical stability:
  # log(sum(w_j * exp(-0.5 * q2_j))) = m + log(sum(w_j * exp(-0.5 * q2_j - m)))
  # where m = max(-0.5 * q2)
  x      <- -0.5 * q2
  shift  <- max(x)
  log_sum_exp_w <- shift + log(sum(weights * exp(x - shift)))

  0.5 * sum(q1) + n1 * log_sum_exp_w
}
