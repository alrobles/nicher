#' Multivariate Gaussian Kernel Density Estimator (Base R Implementation)
#'
#' Computes a multivariate Gaussian kernel density estimate evaluated at a set
#' of points \code{x}, using a sample of reference points \code{data}. This
#' implementation uses only base R, making it fully portable with no external
#' dependencies. The kernel is multivariate normal with covariance matrix
#' \code{H}. If \code{H} is not provided, a diagonal bandwidth matrix based on
#' Scott's rule-of-thumb is used.
#'
#' @param x A numeric matrix of size \code{n_eval x p}, where each row is a point
#'   in \eqn{R^p} at which the density should be evaluated.
#' @param data A numeric matrix of size \code{n_data x p} containing the sample
#'   used to estimate the density.
#' @param H Optional numeric \code{p x p} bandwidth matrix. If \code{NULL}
#'   (default), a diagonal bandwidth matrix is computed using Scott's rule.
#'
#' @return A numeric vector of length \code{n_eval} containing the estimated
#'   density values at each row of \code{x}.
#'
#' @details
#' This function implements a standard multivariate Gaussian kernel estimator:
#' \deqn{
#'   \hat{f}(x) = \frac{1}{n} \sum_{i=1}^n
#'   (2\pi)^{-p/2} |H|^{-1/2}
#'   \exp\left( -\frac{1}{2}(x - x_i)^T H^{-1} (x - x_i) \right).
#' }
#' The default bandwidth matrix is:
#' \deqn{
#'   H = \mathrm{diag}( s^2 \, n^{-2/(p+4)} ),
#' }
#' where \eqn{s} is the sample standard deviation per dimension.
#'
#' @note
#' This implementation is \eqn{O(n_{\text{eval}} \cdot n_{\text{data}})}, and is
#' intended for moderate-sized datasets. For very large problems, a C++ version
#' or tree-based approximation may be preferable.
#'
#' @examples
#' set.seed(1)
#' data <- matrix(rnorm(200), ncol = 2)
#' x    <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
#' kde_gaussian(x, data)
#'
#' @export
kde_gaussian <- function(x, data, H = NULL) {
  x    <- as.matrix(x)
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  
  # Bandwidth por defecto (Scott) en diagonal
  if (is.null(H)) {
    s  <- apply(data, 2, stats::sd)
    bw <- s * n^(-1/(p + 4))
    H  <- diag(bw^2, p)
  }
  
  Hinv <- solve(H)
  detH <- det(H)
  const <- (2 * pi)^(-p/2) * detH^(-1/2)
  
  # Evaluación punto a punto (simple y claro)
  dens <- numeric(nrow(x))
  for (i in seq_len(nrow(x))) {
    dx <- t(t(data) - x[i, ])
    q  <- rowSums((dx %*% Hinv) * dx)
    dens[i] <- const * mean(exp(-0.5 * q))
  }
  dens
}