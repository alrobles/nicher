# Build Cholesky factor L directly from C-vine partial correlations (LKJ).
# - Uses eta > 0 to set per-level Beta shapes: phi_k = eta + (d - k - 1)/2
# - v supplies d*(d-1)/2 reals, one per C-vine edge (in level-major order)
# - Returns lower-triangular L with diag(L) > 0 such that R = L %*% t(L)
#
# Reference: Lewandowski, Kurowicka & Joe (2009), Sec. 2.4 (C-vine) and Eq. (2).

cvine_cholesky <- function(v, d, eta = 1) {
  stopifnot(
    "d must be >= 1" = d >= 1,
    "eta must be > 0" = eta > 0,
    "v must be numeric" = is.numeric(v)
  )

  if (d == 1L) return(matrix(1, 1, 1))

  # Need exactly one real per edge: m = d*(d-1)/2
  m_needed <- d * (d - 1L) / 2L
  if (length(v) != m_needed) {
    stop(sprintf("v must have length %d for d=%d (got %d).",
                 m_needed, d, length(v)))
  }

  cvine_cholesky_cpp(v, d, eta)
}
