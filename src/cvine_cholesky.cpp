#include <Rcpp.h>
using namespace Rcpp;

#include <cmath>
#include <algorithm>

// Forward substitution: solve L_sub * g = b, where L_sub is the top-left
// (j x j) submatrix of the lower-triangular matrix L.
// Used to place row j of the Cholesky factor from the correlation vector r.
static NumericVector forward_solve_lower(const NumericMatrix& L,
                                         const NumericVector& b,
                                         int j) {
  NumericVector g(j);
  for (int i = 0; i < j; i++) {
    double s = b[i];
    for (int k = 0; k < i; k++) {
      s -= L(i, k) * g[k];
    }
    // L(i, i) > 0 is guaranteed because we initialise as identity and set
    // diagonal entries to sqrt(1 - ||g||^2) which stays positive for a
    // positive-definite correlation matrix.
    g[i] = s / L(i, i);
  }
  return g;
}

// cvine_cholesky
//
// Build the Cholesky factor L of a correlation matrix directly from C-vine
// partial correlations following Lewandowski, Kurowicka & Joe (2009), Sec 2.4.
//
// Algorithm overview
// ------------------
// 1. Map each unconstrained real v[idx] -> unit interval via the logistic
//    (sigmoid) function, then to a partial correlation on (-1, 1) via the
//    symmetric Beta quantile:
//      p[k, ell] = 2 * qbeta(sigmoid(v[idx]), phi_k, phi_k) - 1
//    where phi_k = eta + (d - k - 1) / 2  (k in 1 .. d-1, 0-indexed here).
//
// 2. Convert the table of partial correlations to unconditional correlations
//    using the Yule–Kendall (vine) recursion:
//      rho_{ij | L\{m}} = rho_{ij|L} * sqrt((1-r_im^2)(1-r_jm^2)) + r_im*r_jm
//
// 3. For each new row j, solve L[0:j-1, 0:j-1] * g = r (forward substitution)
//    and set L[j, j] = sqrt(1 - ||g||^2).
//
// Parameters
// ----------
// v   : numeric vector of length d*(d-1)/2 (one unconstrained real per C-vine
//       edge, in level-major order, i.e., level k then edges within level k).
// d   : dimension of the target correlation matrix.
// eta : LKJ shape parameter (eta > 0). eta = 1 gives the uniform distribution
//       over correlation matrices; larger eta concentrates mass near identity.
//
// Returns
// -------
// A d x d lower-triangular matrix L such that R = L %*% t(L) is a valid
// correlation matrix with det(R) proportional to det(R)^(eta-1).
//
// Reference: Lewandowski D, Kurowicka D, Joe H (2009). "Generating random
// correlation matrices based on vines and extended onion method." Journal of
// Multivariate Analysis, 100(9):1989-2001.

// [[Rcpp::export(cvine_cholesky_cpp)]]
NumericMatrix cvine_cholesky(NumericVector v, int d, double eta = 1.0) {

  // --- Input validation -----------------------------------------------------
  if (d < 1)   Rcpp::stop("d must be >= 1");
  if (eta <= 0) Rcpp::stop("eta must be > 0 (got %g)", eta);

  // Trivial 1 x 1 case
  if (d == 1) {
    NumericMatrix L(1, 1);
    L(0, 0) = 1.0;
    return L;
  }

  int m_needed = d * (d - 1) / 2;
  if (v.size() != m_needed) {
    Rcpp::stop("v must have length %d for d=%d (got %d).",
               m_needed, d, (int)v.size());
  }

  // --- Step 1: fill partial-correlation table p ---------------------------
  // p(k, ell) = rho_{k+1, ell+1 | 1:k}  for 0 <= k < ell <= d-1  (0-indexed).
  // We use the symmetric Beta distribution on (-1, 1) with shape phi_k.
  // The logistic function provides a stable map from the reals to (0, 1).

  const double eps  = 1e-12;   // clamping for the unit interval
  const double eps15 = 1e-15;  // clamping for the open interval (-1, 1)

  NumericMatrix p(d, d); // initialised to 0

  int idx = 0; // linear index into v
  for (int k = 0; k < d - 1; k++) {
    // phi_k corresponds to R's phi_k = eta + (d - k - 1) / 2  (1-indexed k
    // corresponds to 0-indexed k+1, so substituting: eta + (d-(k+1)-1)/2
    //                                                = eta + (d-k-2)/2 ).
    double phi_k = eta + (d - k - 2) / 2.0;

    for (int ell = k + 1; ell < d; ell++) {
      // Stable sigmoid: u in (eps, 1-eps)
      double u = 1.0 / (1.0 + std::exp(-v[idx++]));
      u = std::min(std::max(u, eps), 1.0 - eps);

      // Map to (-1, 1) via symmetric Beta quantile, then clamp
      double val = 2.0 * R::qbeta(u, phi_k, phi_k, 1, 0) - 1.0;
      p(k, ell) = std::min(std::max(val, -1.0 + eps15), 1.0 - eps15);
    }
  }

  // --- Step 2 & 3: build L row by row ---------------------------------------
  // Initialise L as the d x d identity matrix.
  NumericMatrix L(d, d);
  for (int i = 0; i < d; i++) L(i, i) = 1.0;

  // j (0-indexed) corresponds to R's j = j+1, running from 2..d.
  for (int j = 1; j < d; j++) {

    // Compute the unconditional-correlation vector r of variable j+1 with
    // all previous variables 1..j (0-indexed: 0..j-1).
    NumericVector r(j);

    // r[0] = rho_{1, j+1} = p[0, j]  (direct partial at the first level)
    r[0] = p(0, j);

    // For i = 1..j-1 (R's i = 2..(j-1)), "peel back" conditioning via the
    // Yule–Kendall recursion to recover unconditional correlations.
    for (int i = 1; i < j; i++) {
      // Start from the highest-order partial: rho_{i+1, j+1 | 1:i} = p[i, j]
      double rij = p(i, j);

      // Remove conditioning variables m = i-1, i-2, ..., 0 one by one.
      // Each step converts rho_{i+1,j+1 | {1..i}\{m+1}} using:
      //   rho_{ab|L\{c}} = rho_{ab|L} * sqrt((1-r_ac^2)(1-r_bc^2)) + r_ac*r_bc
      // where r_ac = p[m, i+1] = p(m, i) and r_bc = p[m, j+1] = p(m, j).
      for (int m = i - 1; m >= 0; m--) {
        double rim = p(m, i);
        double rjm = p(m, j);
        double inner = (1.0 - rim * rim) * (1.0 - rjm * rjm);
        rij = rij * std::sqrt(std::max(0.0, inner)) + rim * rjm;
      }
      r[i] = rij;
    }

    // Solve L[0:j-1, 0:j-1] * g = r using forward substitution.
    // g gives the first j entries of the new row j.
    NumericVector g = forward_solve_lower(L, r, j);

    // Fill row j of L
    for (int col = 0; col < j; col++) {
      L(j, col) = g[col];
    }

    // Diagonal entry: maintains ||row j||^2 = 1  (unit-norm row constraint)
    double g_sq_sum = 0.0;
    for (int col = 0; col < j; col++) g_sq_sum += g[col] * g[col];
    L(j, j) = std::sqrt(std::max(0.0, 1.0 - g_sq_sum));
  }

  return L;
}
