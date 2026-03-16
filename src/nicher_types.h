#ifndef NICHER_TYPES_H
#define NICHER_TYPES_H

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

namespace nicher {

// KDE functions
Eigen::VectorXd kde_eigen(const Eigen::MatrixXd& x, const Eigen::MatrixXd& data);
Eigen::VectorXd kde_2d(const Eigen::MatrixXd& x, const Eigen::MatrixXd& data);

// Mahalanobis functions
double sum_mahalanobis_sq(const Eigen::MatrixXd& X,
                          const Eigen::VectorXd& mu,
                          const Eigen::MatrixXd& L);
Eigen::VectorXd mahalanobis_sq_vec(const Eigen::MatrixXd& X,
                                   const Eigen::VectorXd& mu,
                                   const Eigen::MatrixXd& L);

} // namespace nicher

#endif