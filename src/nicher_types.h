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
double loglik_niche_math_kde_bias_corrected_eigen(
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
double loglik_niche_math_kde_bias_corrected_grad_eigen(
    const double* theta, int n_theta,
    const Eigen::MatrixXd& env_occ,
    const Eigen::MatrixXd& M_den,
    const Eigen::VectorXd& w_occ,
    const Eigen::VectorXd& w_den,
    double eta,
    double gradstep_rel, double gradstep_abs,
    double* g_out);

// Math-scale penalized weighted negative log-likelihood. Implements the
// Jimenez & Soberon 2022 (Ecological Modelling 438:109982) Eq. 5 formula
// EXACTLY, plus a weakly-informative ridge penalty on log_sigma that
// stabilises the optimization against the well-known Patil & Ord (1976)
// sigma -> infinity drift in the pure-ML weighted-distribution model.
//
//   -log L = 0.5 sum_i q1(x_i)
//          - sum_i log w(x_i)
//          + n_occ * log( sum_j w(y_j) * exp(-q2(y_j)/2) )
//          + lambda * sum_k (log_sigma_k - log_sigma_center_k)^2
//
// `prior_log_sigma_center` has length p; `prior_log_sigma_lambda` is a
// non-negative scalar (lambda = 0 reduces to plain paper Eq. 5).
double loglik_niche_math_weighted_eigen(
    const double* theta, int n_theta,
    const Eigen::MatrixXd& env_occ,
    const Eigen::MatrixXd& M_den,
    const Eigen::VectorXd& w_occ,
    const Eigen::VectorXd& w_den,
    double eta,
    const Eigen::VectorXd& prior_log_sigma_center,
    double prior_log_sigma_lambda);

// Hybrid analytic gradient for the penalized weighted kernel; same hybrid
// strategy as loglik_niche_math_kde_bias_corrected_grad_eigen plus the closed-form
// gradient of the ridge term on log_sigma.
double loglik_niche_math_weighted_grad_eigen(
    const double* theta, int n_theta,
    const Eigen::MatrixXd& env_occ,
    const Eigen::MatrixXd& M_den,
    const Eigen::VectorXd& w_occ,
    const Eigen::VectorXd& w_den,
    double eta,
    const Eigen::VectorXd& prior_log_sigma_center,
    double prior_log_sigma_lambda,
    double gradstep_rel, double gradstep_abs,
    double* g_out);

// Math-scale skew-normal presence-only negative log-likelihood. Theta layout:
//   [mu(0..p-1), log_sigma(0..p-1), v(0..p*(p-1)/2 - 1), alpha(0..p-1)]
// (the additional alpha block extends the Gaussian theta by p parameters).
//
// Density: f_{SN}(x; mu, Sigma, alpha) = 2 * phi_k(x-mu; Sigma) * Phi(z(x))
// where z(x) = sum_k alpha_k * (x_k - mu_k) / sigma_k.
// Reference: Azzalini & Capitanio (1999), J. R. Stat. Soc. Ser. B 61(3).
double loglik_niche_math_skew_normal_eigen(
    const double* theta, int n_theta,
    const Eigen::MatrixXd& env_occ, double eta);

// Math-scale skew-normal weighted negative log-likelihood. Theta layout
// matches loglik_niche_math_skew_normal_eigen above. Composes the
// skew-normal base with paper Eq. 5 of Jimenez & Soberon (2022) plus a
// ridge penalty on log_sigma:
//
//   -log L = 0.5 sum_i q1(x_i) - sum_i log Phi(z1_i)
//          - sum_i log w(x_i)
//          + n_occ * log( sum_j w(y_j) * exp(-q2(y_j)/2) * Phi(z2_j) )
//          + lambda * sum_k (log_sigma_k - log_sigma_center_k)^2
//
// (Constants log 2 and (2 pi)^{-k/2} |Sigma|^{-1/2} cancel between
// numerator and denominator log-sum-exp, just as in the Gaussian case.)
double loglik_niche_math_skew_normal_weighted_eigen(
    const double* theta, int n_theta,
    const Eigen::MatrixXd& env_occ,
    const Eigen::MatrixXd& M_den,
    const Eigen::VectorXd& w_occ,
    const Eigen::VectorXd& w_den,
    double eta,
    const Eigen::VectorXd& prior_log_sigma_center,
    double prior_log_sigma_lambda);

// Math-scale skew-t presence-only (multivariate non-central skew-t via 1-D
// Gauss-Laguerre quadrature on the chi^2_r mixing variable). Theta layout:
//   [mu(0..p-1), log_sigma(0..p-1), v(0..p*(p-1)/2 - 1),
//    alpha(0..p-1), log_r(0)]
// (extends the skew-normal layout by 1 extra parameter at the end: log r).
//
// Reference: Branco & Dey (2001), J. Multivariate Analysis 79(1):99-113.
// Density:
//   f_T(t) = 2 (2*pi)^{-k/2} |Sigma|^{-1/2} r^{-k/2} / [2^{r/2} Gamma(r/2)]
//          * (2 / A(t))^{(k+r)/2} * I(t)
//   A(t)   = 1 + q(t)/r,   q(t) = (t-mu)^T Sigma^{-1} (t-mu)
//   I(t)   = integral_0^inf u^{(k+r)/2 - 1} exp(-u)
//                          * Phi(z(t) sqrt(2 u / (A(t) r))) du
// approximated by Q = 32 standard Gauss-Laguerre nodes/weights.
//
// As r -> infinity the kernel reduces to skew_normal modulo a constant
// in theta (verified analytically and numerically in test-loglik_math_cpp_parity).
double loglik_niche_math_skew_t_eigen(
    const double* theta, int n_theta,
    const Eigen::MatrixXd& env_occ, double eta);

// Math-scale skew-t weighted (NCST density + paper Eq. 5 + ridge prior on
// log_sigma). Same theta layout as loglik_niche_math_skew_t_eigen.
//
// Constants C(r, Sigma) = 2 (2*pi)^{-k/2} |Sigma|^{-1/2} r^{-k/2} /
//                         [2^{r/2} Gamma(r/2)] cancel exactly between
// numerator and denominator log-sum-exp (just like the skew-normal case),
// so the weighted form does not involve log|Sigma|, log Gamma(r/2), or
// log(r) prefactors directly:
//
//   -log L = sum_i ((k+r)/2) log A(t_i) - sum_i log I(t_i)
//          - sum_i log w(t_i)
//          + n_occ * log( sum_j w(y_j) A(y_j)^{-(k+r)/2} I(y_j) )
//          + lambda * sum_k (log_sigma_k - log_sigma_center_k)^2
double loglik_niche_math_skew_t_weighted_eigen(
    const double* theta, int n_theta,
    const Eigen::MatrixXd& env_occ,
    const Eigen::MatrixXd& M_den,
    const Eigen::VectorXd& w_occ,
    const Eigen::VectorXd& w_den,
    double eta,
    const Eigen::VectorXd& prior_log_sigma_center,
    double prior_log_sigma_lambda);

}

#endif