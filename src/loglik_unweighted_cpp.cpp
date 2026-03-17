// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Log-likelihood for the unweighted Gaussian niche model (C++ engine)
//'
//' Computes the Gaussian log-likelihood for the unweighted niche model using
//' Armadillo linear algebra routines.  Equivalent to \code{loglik_unweighted_math}
//' but implemented in C++ for performance.
//'
//' @param sam1 Numeric matrix of presence points (rows = observations, columns = variables).
//' @param sam2 Numeric matrix of background environmental points.
//' @param mu Numeric vector of means (optimum) of length p.
//' @param L Lower triangular Cholesky factor matrix (p x p). Diagonal entries must be positive.
//' @return Log-likelihood value (scalar double).
//' @export
// [[Rcpp::export]]
double loglik_unweighted_cpp(
    const arma::mat& sam1,
    const arma::mat& sam2,
    const arma::vec& mu,
    const arma::mat& L
) {
  arma::mat S = L * L.t();

  // log|S| = 2 * sum(log(diag(L)))
  double logdet = 2.0 * arma::sum(arma::log(arma::diagvec(L)));

  arma::mat S_inv = arma::inv_sympd(S);

  int n1 = sam1.n_rows;
  int n2 = sam2.n_rows;

  double sum_q1 = 0.0;
  for (int i = 0; i < n1; i++) {
    arma::vec d = sam1.row(i).t() - mu;
    sum_q1 += arma::as_scalar(d.t() * S_inv * d);
  }

  double sum_q2 = 0.0;
  for (int i = 0; i < n2; i++) {
    arma::vec d = sam2.row(i).t() - mu;
    sum_q2 += arma::as_scalar(d.t() * S_inv * d);
  }

  double logL = -0.5 * (n1 * logdet + sum_q1) + 0.5 * (n2 * logdet + sum_q2);
  // Note: the sign difference is intentional — the log-likelihood is computed as
  // the difference of average log-densities (presence minus background), which
  // gives a positive contribution from the background term.
  return logL;
}
