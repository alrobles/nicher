#' Multivariate Gaussian KDE with fixed Scott bandwidth
#'
#' Computes a multivariate Gaussian kernel density estimate at the rows of `x`,
#' using `data` as the reference sample. Bandwidths are set internally via
#' Scott's rule-of-thumb (diagonal). For 2D, a manually optimized version is used;
#' for higher dimensions, an Eigen-based vectorized version is used.
#'
#' @param x Numeric matrix `n_eval x p` (evaluation points).
#' @param data Numeric matrix `n_data x p` (reference sample for KDE).
#'
#' @return Numeric vector of length `n_eval` with the density estimates.
#' @examples
#' x    <- as.matrix(example_env_occ_2d)
#' data <- as.matrix(example_env_m_2d)
#' dens <- kde_gaussian(x, data)
#' head(dens)
#' @export
kde_gaussian <- function(x, data) {
  x <- as.matrix(x)
  data <- as.matrix(data)
  if (ncol(x) != ncol(data)) {
    stop("x and data must have the same number of columns")
  }
  p <- ncol(x)
  if (p == 2) {
    kde_gaussian_2d_cpp(x, data)
  } else {
    kde_gaussian_eigen_cpp(x, data)
  }
}
