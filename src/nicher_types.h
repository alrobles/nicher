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

// Eigen-native C-vine Cholesky factor of a correlation matrix
// (Lewandowski-Kurowicka-Joe 2009). Writes a lower-triangular d x d matrix
// L_out such that R = L_out * L_out^T is a valid correlation matrix. The
// vector v has length d*(d-1)/2 of unconstrained reals (level-major).
void cvine_cholesky_eigen(const Eigen::Ref<const Eigen::VectorXd>& v,
                          int d, double eta,
                          Eigen::MatrixXd& L_out);

// Math-scale presence-only negative log-likelihood. Theta layout:
//   [mu(0..p-1), log_sigma(0..p-1), v(0..p*(p-1)/2 - 1)].
double loglik_niche_math_presence_only_eigen(
    const double* theta, int n_theta,
    const Eigen::MatrixXd& env_occ, double eta);

// Math-scale weighted negative log-likelihood with PRECOMPUTED KDE weights.
// w_occ[i] is the KDE weight at presence row i; w_den[j] at denominator row j.
double loglik_niche_math_weighted_eigen(
    const double* theta, int n_theta,
    const Eigen::MatrixXd& env_occ,
    const Eigen::MatrixXd& M_den,
    const Eigen::VectorXd& w_occ,
    const Eigen::VectorXd& w_den,
    double eta);

// Hybrid analytic gradient for the weighted kernel:
//   - mu and log_sigma blocks: closed-form
//   - v block (C-vine partials): central finite difference
// Writes f and the n_theta-length gradient vector g_out. Returns f.
double loglik_niche_math_weighted_grad_eigen(
    const double* theta, int n_theta,
    const Eigen::MatrixXd& env_occ,
    const Eigen::MatrixXd& M_den,
    const Eigen::VectorXd& w_occ,
    const Eigen::VectorXd& w_den,
    double eta,
    double gradstep_rel, double gradstep_abs,
    double* g_out);

}

#endif