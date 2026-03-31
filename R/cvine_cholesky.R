#' Build Cholesky factor of a correlation matrix from C‑vine partial correlations
#'
#' Constructs the lower‑triangular Cholesky factor \eqn{\mathbf{L}} of a
#' \eqn{d \times d} correlation matrix using the C‑vine method described in
#' Lewandowski, Kurowicka & Joe (2009). The factor satisfies
#' \eqn{\mathbf{R} = \mathbf{L}\mathbf{L}^\top} with \eqn{\operatorname{diag}(\mathbf{R}) = 1}.
#'
#' @param v Numeric vector of unconstrained reals, one for each C‑vine edge.
#'   The length must be \eqn{d(d-1)/2}. For \eqn{d = 2}, a single value is required.
#'   The order follows a **level‑major** sequence: first all edges of level 1
#'   (i.e., between variable 1 and each later variable), then edges of level 2
#'   (between variable 2 and later variables, conditioned on variable 1), and so on.
#' @param d Integer, dimension of the target correlation matrix (\eqn{d \ge 1}).
#' @param eta Positive numeric shape parameter for the LKJ‑C‑vine prior
#'   (default \code{1}). \eqn{\eta = 1} gives a uniform distribution over
#'   correlation matrices; larger values concentrate mass near the identity.
#'
#' @return A \eqn{d \times d} lower‑triangular matrix \eqn{\mathbf{L}} with
#'   positive diagonal entries such that \eqn{\mathbf{L}\mathbf{L}^\top} is a
#'   valid correlation matrix (unit diagonal, positive definite). For \eqn{d = 1},
#'   returns a \eqn{1 \times 1} matrix with entry \code{1}.
#'
#' @details
#' The algorithm proceeds in three steps:
#' \enumerate{
#'   \item Each element of \code{v} is mapped to the unit interval via the
#'         logistic (sigmoid) function, then transformed to a partial correlation
#'         on \eqn{(-1, 1)} using the quantile function of a symmetric Beta
#'         distribution with shape \eqn{\phi_k = \eta + (d - k - 2)/2}
#'         (for level \eqn{k}, 0‑indexed).
#'   \item The table of partial correlations is converted to unconditional
#'         correlations using the Yule–Kendall recursion (vine recursion).
#'   \item For each new row \eqn{j} (starting from \eqn{j = 2}), the algorithm
#'         solves a triangular system to obtain the first \eqn{j-1} entries of the
#'         row, and sets the diagonal entry to maintain unit row norm.
#' }
#' The resulting \eqn{\mathbf{L}} can be used directly to construct the correlation
#' matrix (\code{tcrossprod(L)}), or to build a covariance matrix by scaling
#' with standard deviations.
#'
#' @references
#' Lewandowski, D., Kurowicka, D., & Joe, H. (2009).
#' Generating random correlation matrices based on vines and extended onion method.
#' *Journal of Multivariate Analysis*, 100(9), 1989–2001.
#' \doi{10.1016/j.jmva.2009.04.008}
#'
#'
#' @examples
#' # For a 2x2 correlation matrix, we need one parameter
#' v <- 0.5
#' L <- cvine_cholesky(v, d = 2, eta = 1)
#' R <- tcrossprod(L)
#' print(R)
#'
#' # For 3 dimensions, we need 3 parameters (d*(d-1)/2 = 3)
#' v <- c(0.1, -0.2, 0.8)
#' L <- cvine_cholesky(v, d = 3)
#' R <- tcrossprod(L)
#' all.equal(diag(R), rep(1, 3)) # Should be TRUE
#'
#' @export
cvine_cholesky <- function(v, d, eta = 1) {
  stopifnot(
    "d must be >= 1" = d >= 1,
    "eta must be > 0" = eta > 0,
    "v must be numeric" = is.numeric(v)
  )

  if (d == 1L) {
    return(matrix(1, 1, 1))
  }

  # Need exactly one real per edge: m = d*(d-1)/2
  m_needed <- d * (d - 1L) / 2L
  if (length(v) != m_needed) {
    stop(sprintf(
      "v must have length %d for d=%d (got %d).",
      m_needed, d, length(v)
    ))
  }

  cvine_cholesky_cpp(v, d, eta)
}
