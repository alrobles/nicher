#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector kde_gaussian_eigen_cpp(const NumericMatrix& x, const NumericMatrix& data) {
  int n_eval = x.nrow();
  int n_data = data.nrow();
  int p = data.ncol();
  
  if (n_data <= 0 || p <= 0) stop("data must have n>0 and p>0");
  if (x.ncol() != p) stop("x and data must have the same number of columns");
  
  Eigen::Map<Eigen::MatrixXd> X(as<Eigen::Map<Eigen::MatrixXd> >(x));
  Eigen::Map<Eigen::MatrixXd> D(as<Eigen::Map<Eigen::MatrixXd> >(data));
  
  // Scott's rule: h_j = sd_j * n^(-1/(p+4))
  Eigen::VectorXd s(p);
  for (int j = 0; j < p; ++j) {
    Eigen::VectorXd col = D.col(j);
    double mean = col.mean();
    double var = (col.array() - mean).square().mean();
    s(j) = std::sqrt(var);
    if (s(j) <= 0) s(j) = 1e-8;
  }
  
  double nf = std::pow(static_cast<double>(n_data), -1.0 / (p + 4.0));
  Eigen::VectorXd bw = s * nf;
  Eigen::VectorXd inv_bw2 = (bw.array() * bw.array()).inverse();
  double prod_bw = bw.prod();
  if (prod_bw <= 0) prod_bw = 1e-32;
  
  double log2pi = std::log(2.0 * M_PI);
  double log_const = -0.5 * p * log2pi - std::log(prod_bw);
  double const_kernel = std::exp(log_const);
  
  // Escalar matrices
  Eigen::MatrixXd X_scaled = X * inv_bw2.asDiagonal();
  Eigen::MatrixXd D_scaled = D * inv_bw2.asDiagonal();
  
  Eigen::VectorXd norm_X = X_scaled.rowwise().squaredNorm();
  Eigen::VectorXd norm_D = D_scaled.rowwise().squaredNorm();
  
  Eigen::MatrixXd cross = -2.0 * X_scaled * D_scaled.transpose();
  Eigen::MatrixXd dist2 = cross.colwise() + norm_X;
  dist2.rowwise() += norm_D.transpose();
  
  Eigen::MatrixXd kernel = (-0.5 * dist2).array().exp();
  Eigen::VectorXd dens = const_kernel * (kernel.rowwise().sum() / n_data);
  
  return wrap(dens);
}