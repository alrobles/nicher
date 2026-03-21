#ifndef NICHER_TYPES_H
#define NICHER_TYPES_H

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

namespace nicher {

// Large finite penalty returned when the log-likelihood evaluates to NaN/Inf.
// This allows the optimizer to recover by steering away from the problematic
// region instead of crashing with a type conversion error.
constexpr double OPTIM_PENALTY = 1e300;

// Minimum KDE weight used before taking log() to prevent log(0) = -Inf.
// Set near the smallest representable positive double to avoid distorting
// the likelihood surface while preventing -Inf/NaN propagation.
constexpr double MIN_KDE_WEIGHT = 1e-300;

Eigen::VectorXd kde_2d(const Eigen::MatrixXd& x, const Eigen::MatrixXd& data);
Eigen::VectorXd kde_eigen(const Eigen::MatrixXd& x, const Eigen::MatrixXd& data);
double sum_mahalanobis_sq(const Eigen::MatrixXd& X, const Eigen::VectorXd& mu, const Eigen::MatrixXd& L);
Eigen::VectorXd mahalanobis_sq_vec(const Eigen::MatrixXd& X, const Eigen::VectorXd& mu, const Eigen::MatrixXd& L);
}

#endif