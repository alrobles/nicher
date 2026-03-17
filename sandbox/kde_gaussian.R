#' Multivariate Gaussian KDE with fixed Scott bandwidth (C++ backend only)
#'
#' Computes a multivariate Gaussian kernel density estimate at the rows of `x`,
#' using `data` as the reference sample. Bandwidths are set internally via
#' Scott's rule-of-thumb (diagonal), and the computation is performed entirely
#' in C++ (Rcpp).
#'
#' @param x Numeric matrix `n_eval x p` (evaluation points).
#' @param data Numeric matrix `n_data x p` (reference sample for KDE).
#'
#' @return Numeric vector of length `n_eval` with the density estimates.
#' @export
kde_gaussian <- function(x, data) {
  x    <- as.matrix(x)
  data <- as.matrix(data)
  if (ncol(x) != ncol(data)) {
    stop("x and data must have the same number of columns")
  }
  kde_gaussian_rcpp(x, data)
}