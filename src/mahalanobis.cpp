#include "nicher_types.h"

namespace nicher {

double sum_mahalanobis_sq(const Eigen::MatrixXd& X,
                          const Eigen::VectorXd& mu,
                          const Eigen::MatrixXd& L) {
  int n = X.rows();
  int p = X.cols();
  Eigen::MatrixXd diff(p, n);
  for (int j = 0; j < n; ++j) {
    diff.col(j) = X.row(j).transpose() - mu;
  }
  Eigen::MatrixXd y = L.triangularView<Eigen::Lower>().solve(diff);
  return y.colwise().squaredNorm().sum();
}

Eigen::VectorXd mahalanobis_sq_vec(const Eigen::MatrixXd& X,
                                   const Eigen::VectorXd& mu,
                                   const Eigen::MatrixXd& L) {
  int n = X.rows();
  int p = X.cols();
  Eigen::MatrixXd diff(p, n);
  for (int j = 0; j < n; ++j) {
    diff.col(j) = X.row(j).transpose() - mu;
  }
  Eigen::MatrixXd y = L.triangularView<Eigen::Lower>().solve(diff);
  return y.colwise().squaredNorm();
}

} // namespace nicher