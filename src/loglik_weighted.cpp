// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Log-likelihood for the weighted Gaussian niche model (C++ engine)
//'
//' Computes the weighted presence-only (semi-log-likelihood) for a Gaussian
//' niche model using Mahalanobis quadratic forms.  Background points can be
//' given non-negative weights so that denser regions of environmental space
//' contribute proportionally more to the normalizing integral.  Weights are
//' normalized internally so that they sum to one.
//'
//' When all weights are equal the result is identical to
//' \code{loglik_presenceonly_cpp}.
//'
//' @param sam1 Numeric matrix of presence points
//'   (rows = observations, columns = variables).
//' @param sam2 Numeric matrix of background environmental points.
//' @param mu Numeric vector of means (optimum) of length p.
//' @param S Covariance matrix (p x p).
//' @param weights Numeric vector of non-negative weights for the background
//'   points in \code{sam2}.  Length must equal \code{nrow(sam2)}.
//' @return Negative log-likelihood value (scalar double).
//' @export
// [[Rcpp::export]]
double loglik_weighted_cpp(
    const arma::mat& sam1,
    const arma::mat& sam2,
    const arma::vec& mu,
    const arma::mat& S,
    const arma::vec& weights
) {
  if ((int)weights.n_elem != (int)sam2.n_rows) {
    Rcpp::stop("'weights' length must equal nrow(sam2)");
  }
  if (arma::any(weights < 0)) {
    Rcpp::stop("'weights' must be non-negative");
  }
  double w_sum = arma::sum(weights);
  if (w_sum == 0.0) {
    Rcpp::stop("'weights' must not all be zero");
  }
  arma::vec w = weights / w_sum;

  arma::mat S_inv = arma::inv_sympd(S);

  int n1 = sam1.n_rows;
  int n2 = sam2.n_rows;

  double sum_q1 = 0.0;
  for (int i = 0; i < n1; i++) {
    arma::vec d = sam1.row(i).t() - mu;
    sum_q1 += arma::as_scalar(d.t() * S_inv * d);
  }

  arma::vec q2(n2);
  for (int i = 0; i < n2; i++) {
    arma::vec d = sam2.row(i).t() - mu;
    q2(i) = arma::as_scalar(d.t() * S_inv * d);
  }

  // Weighted log-sum-exp for numerical stability:
  // log(sum(w_j * exp(-0.5 * q2_j))) = m + log(sum(w_j * exp(-0.5 * q2_j - m)))
  // where m = max(-0.5 * q2).
  arma::vec x = -0.5 * q2;
  double shift = arma::max(x);
  double log_sum_exp_w = shift + std::log(arma::sum(w % arma::exp(x - shift)));

  return 0.5 * sum_q1 + n1 * log_sum_exp_w;
}
