// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

//' Compute the Silverman bandwidth matrix for multivariate KDE
//'
//' Applies Silverman's rule of thumb to derive a full bandwidth matrix from
//' the sample covariance.  The scalar bandwidth factor is
//' \eqn{h = n^{-2/(p+4)}} and the bandwidth matrix is \eqn{H = h \cdot S},
//' where S is the sample covariance of \code{ref_sample}.
//'
//' This function is intended to be called once on the reference (background)
//' sample so that the result can be cached and reused across many KDE
//' evaluations (e.g., during likelihood optimisation).
//'
//' @param ref_sample Numeric matrix of reference points (rows = observations,
//'   columns = variables).
//' @return Bandwidth matrix H (p x p, symmetric positive definite).
//' @export
//'
//' @examples
//' H <- kde_bandwidth_silverman(as.matrix(samMPts))
// [[Rcpp::export]]
arma::mat kde_bandwidth_silverman(const arma::mat& ref_sample) {
  int n = ref_sample.n_rows;
  int p = ref_sample.n_cols;
  arma::mat S = arma::cov(ref_sample);
  double h = std::pow(static_cast<double>(n), -2.0 / (p + 4));
  return h * S;
}

//' Evaluate multivariate Gaussian KDE at query points (parallelised)
//'
//' Evaluates the log-density of a multivariate Gaussian kernel density
//' estimate at each row of \code{query}.  The expensive inner loop (over
//' reference points) is parallelised with OpenMP when available.
//'
//' Separating bandwidth computation (\code{\link{kde_bandwidth_silverman}})
//' from evaluation enables \emph{caching}: compute \code{H} once for a fixed
//' background sample, then call \code{kde_eval_cpp} repeatedly during
//' optimisation without re-computing the bandwidth.
//'
//' @param query Numeric matrix of evaluation points (m x p).
//' @param ref Numeric matrix of reference (training) points (n x p).
//' @param H Bandwidth matrix (p x p), e.g. from
//'   \code{\link{kde_bandwidth_silverman}}.
//' @return Numeric vector of length m containing log-density values.
//' @export
//'
//' @examples
//' H   <- kde_bandwidth_silverman(as.matrix(samMPts))
//' lkd <- kde_eval_cpp(as.matrix(spOccPnts), as.matrix(samMPts), H)
// [[Rcpp::export]]
arma::vec kde_eval_cpp(
    const arma::mat& query,
    const arma::mat& ref,
    const arma::mat& H
) {
  int m = query.n_rows;
  int n = ref.n_rows;
  int p = ref.n_cols;

  arma::mat H_inv = arma::inv_sympd(H);
  double log_det_H = arma::log_det_sympd(H);
  double log_norm = -0.5 * static_cast<double>(p) *
                    std::log(2.0 * arma::datum::pi) -
                    0.5 * log_det_H;
  double log_n = std::log(static_cast<double>(n));

  arma::vec log_density(m);

  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  // Preallocate one working buffer per thread to avoid per-iteration heap
  // allocation (n doubles per thread rather than n doubles per query point).
  std::vector<std::vector<double>> lw_buf(
      nthreads, std::vector<double>(n));

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < m; i++) {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    std::vector<double>& lw = lw_buf[tid];
    for (int j = 0; j < n; j++) {
      arma::vec d = query.row(i).t() - ref.row(j).t();
      lw[j] = log_norm - 0.5 * arma::as_scalar(d.t() * H_inv * d);
    }
    // Numerically stable log-sum-exp
    double max_lw = *std::max_element(lw.begin(), lw.end());
    double sum_exp = 0.0;
    for (int j = 0; j < n; j++) sum_exp += std::exp(lw[j] - max_lw);
    log_density(i) = max_lw + std::log(sum_exp) - log_n;
  }

  return log_density;
}
