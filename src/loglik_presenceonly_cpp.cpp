// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

//' Log-likelihood for the presence-only niche model (C++ engine)
//'
//' Computes the presence-only (semi-log-likelihood) for a Gaussian niche model
//' using Mahalanobis quadratic forms.  Equivalent to \code{loglik_presenceonly_math}
//' but implemented in C++ for performance.
//'
//' @param sam1 Numeric matrix of presence points (rows = observations, columns = variables).
//' @param sam2 Numeric matrix of background environmental points.
//' @param mu Numeric vector of means of length p.
//' @param S Covariance matrix (p x p).
//' @return Negative log-likelihood value (scalar double).
//' @export
// [[Rcpp::export]]
double loglik_presenceonly_cpp(
    const arma::mat& sam1,
    const arma::mat& sam2,
    const arma::vec& mu,
    const arma::mat& S
) {
  arma::mat S_inv = arma::inv_sympd(S);

  int n1 = sam1.n_rows;
  int n2 = sam2.n_rows;

  double sum_q1 = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:sum_q1) schedule(static)
#endif
  for (int i = 0; i < n1; i++) {
    arma::vec d = sam1.row(i).t() - mu;
    sum_q1 += arma::as_scalar(d.t() * S_inv * d);
  }

  arma::vec q2(n2);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < n2; i++) {
    arma::vec d = sam2.row(i).t() - mu;
    q2(i) = arma::as_scalar(d.t() * S_inv * d);
  }

  // Use log-sum-exp for numerical stability:
  // For x_i = -0.5 * q2(i), log(sum(exp(x_i))) = m + log(sum(exp(x_i - m)))
  // where m = max(x_i) = -0.5 * min(q2).
  double shift = -0.5 * arma::min(q2);
  double log_sum_exp = shift + std::log(arma::sum(arma::exp(-0.5 * q2 - shift)));

  double neg_log = 0.5 * sum_q1 + n1 * log_sum_exp;
  return neg_log;
}
