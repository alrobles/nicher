#' Convert math-scale parameters to biological-scale for the Gaussian niche model
#'
#' @param param_vector Named numeric vector of parameters in mathematical scale.
#'   Names: \code{mu1..mup}, \code{L1..Lk} (k = p*(p+1)/2).
#'   Diagonal entries of L are on the log scale.
#' @return A list with elements \code{mu} (numeric vector) and \code{L} (lower
#'   triangular matrix).
math_to_bio_niche <- function(param_vector) {
  mu <- param_vector[grep("^mu", names(param_vector))]
  p <- length(mu)

  L_elements <- param_vector[grep("^L", names(param_vector))]
  k <- length(L_elements)
  if (k != p * (p + 1) / 2) {
    stop("Incorrect number of L elements")
  }

  L <- matrix(0, p, p)
  L[lower.tri(L, diag = TRUE)] <- L_elements
  diag(L) <- exp(diag(L))

  list(mu = mu, L = L)
}
