#include "nicher_types.h"

// Explicit declaration (optional, aids diagnostics)
namespace nicher {
Eigen::VectorXd kde_2d(const Eigen::MatrixXd& x, const Eigen::MatrixXd& data);
}

// [[Rcpp::export]]
Rcpp::NumericVector kde_gaussian_2d_cpp(const Rcpp::NumericMatrix& x,
                                        const Rcpp::NumericMatrix& data) {
  Eigen::Map<Eigen::MatrixXd> X(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(x));
  Eigen::Map<Eigen::MatrixXd> D(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(data));
  Eigen::VectorXd res = nicher::kde_2d(X, D);
  return Rcpp::wrap(res);
}