#' Negative log likelyhood
#' @param sam1 A data.frame containing the original sample of environmental combinations that correspond to presences
#' @param sam2 A data.frame containing a second random sample of environmental combinations which come from the area of study (M)
#' @param mu A vector mu of parameters
#' @param S The covariance matrix from environmental data frame
#'
#' @return A negative log likelihood value
#' @export
#'
#' @examples
#' par <- get_ellip_par(spOccPnts)
#' negloglike_multivariable(par$mu, par$S, spOccPnts, samMPts)
negloglike_multivariable_matrix <- function(mu, S, sam1, sam2){

  q1 <- stats::mahalanobis(x = sam1, center = mu, cov = S, inverted = FALSE) # quadratic terms of presence points
  q2 <- stats::mahalanobis(x = sam2, center = mu, cov = S, inverted = FALSE) # quadratic terms of M points

  # Check why Cholesky decomposition doesn't works
  # It gets a non symetric matrix error.
  # q1 <- mvnfast::maha(X = sam1, mu = mu, sigma = S, isChol = FALSE) # quadratic terms of presence points
  # q2 <- mvnfast::maha(X = sam2, mu = mu, sigma = S, isChol = FALSE) # quadratic terms of M points


  # negative log-likelihood value
  neg_log <- get_negative_log(q1, q2)
  return(neg_log)
}
