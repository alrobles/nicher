#include "nicher_types.h"

// [[Rcpp::export]]
Rcpp::NumericVector kde_gaussian_eigen_cpp(const Rcpp::NumericMatrix& x,
                                           const Rcpp::NumericMatrix& data) {
  Eigen::Map<Eigen::MatrixXd> X(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(x));
  Eigen::Map<Eigen::MatrixXd> D(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(data));
  Eigen::VectorXd res = nicher::kde_eigen(X, D);
  return Rcpp::wrap(res);
}