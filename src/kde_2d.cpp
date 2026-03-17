#include "nicher_types.h"

namespace nicher {

Eigen::VectorXd kde_2d(const Eigen::MatrixXd& x, const Eigen::MatrixXd& data) {
  int n_eval = x.rows();
  int n_data = data.rows();
  
  double mean1 = data.col(0).mean();
  double mean2 = data.col(1).mean();
  
  double var1 = (data.col(0).array() - mean1).square().mean();
  double var2 = (data.col(1).array() - mean2).square().mean();
  
  double s1 = std::sqrt(var1);
  double s2 = std::sqrt(var2);
  if (s1 <= 0) s1 = 1e-8;
  if (s2 <= 0) s2 = 1e-8;
  
  double nf = std::pow(static_cast<double>(n_data), -1.0 / 6.0);
  double bw1 = s1 * nf;
  double bw2 = s2 * nf;
  double inv_bw1_2 = 1.0 / (bw1 * bw1);
  double inv_bw2_2 = 1.0 / (bw2 * bw2);
  double prod_bw = bw1 * bw2;
  if (prod_bw <= 0) prod_bw = 1e-32;
  
  double log2pi = std::log(2.0 * M_PI);
  double log_const = -0.5 * 2 * log2pi - std::log(prod_bw);
  double const_kernel = std::exp(log_const);
  
  Eigen::VectorXd dens(n_eval);
  for (int i = 0; i < n_eval; ++i) {
    double xi1 = x(i, 0);
    double xi2 = x(i, 1);
    double sum = 0.0;
    for (int j = 0; j < n_data; ++j) {
      double d1 = data(j, 0) - xi1;
      double d2 = data(j, 1) - xi2;
      double q = d1 * d1 * inv_bw1_2 + d2 * d2 * inv_bw2_2;
      sum += std::exp(-0.5 * q);
    }
    dens(i) = const_kernel * (sum / n_data);
  }
  return dens;
}

} // namespace nicher