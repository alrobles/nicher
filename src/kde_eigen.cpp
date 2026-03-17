#include "nicher_types.h"

namespace nicher {

Eigen::VectorXd kde_eigen(const Eigen::MatrixXd& x, const Eigen::MatrixXd& data) {
  int n_eval = x.rows();
  int n_data = data.rows();
  int p = data.cols();
  
  // Scott's rule
  Eigen::VectorXd s(p);
  for (int j = 0; j < p; ++j) {
    Eigen::VectorXd col = data.col(j);
    double mean = col.mean();
    double var = (col.array() - mean).square().mean();
    s(j) = std::sqrt(var);
    if (s(j) <= 0) s(j) = 1e-8;
  }
  
  double nf = std::pow(static_cast<double>(n_data), -1.0 / (p + 4.0));
  Eigen::VectorXd bw = s * nf;
  Eigen::VectorXd inv_bw = bw.array().inverse();
  double prod_bw = bw.prod();
  if (prod_bw <= 0) prod_bw = 1e-32;
  
  double log2pi = std::log(2.0 * M_PI);
  double log_const = -0.5 * p * log2pi - std::log(prod_bw);
  double const_kernel = std::exp(log_const);
  
  Eigen::MatrixXd X_scaled = x * inv_bw.asDiagonal();
  Eigen::MatrixXd D_scaled = data * inv_bw.asDiagonal();
  
  Eigen::VectorXd norm_X = X_scaled.rowwise().squaredNorm();
  Eigen::VectorXd norm_D = D_scaled.rowwise().squaredNorm();
  
  Eigen::MatrixXd cross = -2.0 * X_scaled * D_scaled.transpose();
  Eigen::MatrixXd dist2 = cross.colwise() + norm_X;
  dist2.rowwise() += norm_D.transpose();
  
  Eigen::MatrixXd kernel = (-0.5 * dist2).array().exp();
  Eigen::VectorXd dens = const_kernel * (kernel.rowwise().sum() / n_data);
  return dens;
}

} // namespace nicher