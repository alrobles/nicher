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

// ---------------------------------------------------------------------------
// PriorParams: shared penalised-MLE configuration.
//
// Encapsulates the three optional ridge penalties applied on top of the
// base negative log-likelihood:
//
//   penalty(theta) =
//     mu_lambda        * sum_k ((mu_k - mu_center_k) / exp(log_sigma_center_k))^2
//   + log_sigma_lambda * sum_k (log_sigma_k - log_sigma_center_k)^2
//   + alpha_lambda     * sum_k alpha_k^2                  (skew families only)
//
// The mu penalty is divided by exp(log_sigma_center_k) so that mu_lambda is
// unit-free across variables whose scales differ wildly (e.g. bio1 in degC
// vs bio12 in mm). With the default anchor at the presence-only fit,
// mu_lambda = 1 means "mu is allowed to drift ~1 standard deviation from
// mu_PO before the penalty pushes back".
//
// All lambdas default to 0 (no penalty); when every lambda is 0 the kernel
// reduces exactly to its un-penalised form.
//
// alpha_lambda is ignored if the active kernel does not have an alpha block.
// alpha is always shrunk toward 0 (alpha = 0 recovers the Gaussian sub-model).
// ---------------------------------------------------------------------------
struct PriorParams {
  Eigen::VectorXd mu_center;          // length p (zeros if mu_lambda = 0)
  double          mu_lambda        = 0.0;
  Eigen::VectorXd log_sigma_center;   // length p (zeros if both sigma- and
                                      // mu-lambdas are 0)
  double          log_sigma_lambda = 0.0;
  double          alpha_lambda     = 0.0;
};

// Compute the additive penalty value for a given theta. theta layout:
//   [mu(p), log_sigma(p), v(n_v), alpha(p)?, log_r(1)?]
// alpha block is only present if has_alpha == true.
inline double prior_penalty_value(const double* theta, int p, int n_v,
                                  bool has_alpha,
                                  const PriorParams& pp) {
  double pen = 0.0;
  if (pp.mu_lambda > 0.0 && pp.mu_center.size() == p &&
      pp.log_sigma_center.size() == p) {
    double s = 0.0;
    for (int k = 0; k < p; ++k) {
      const double scale = std::exp(pp.log_sigma_center(k));
      const double d = (theta[k] - pp.mu_center(k)) / scale;
      s += d * d;
    }
    pen += pp.mu_lambda * s;
  }
  if (pp.log_sigma_lambda > 0.0 && pp.log_sigma_center.size() == p) {
    double s = 0.0;
    for (int k = 0; k < p; ++k) {
      const double d = theta[p + k] - pp.log_sigma_center(k);
      s += d * d;
    }
    pen += pp.log_sigma_lambda * s;
  }
  if (pp.alpha_lambda > 0.0 && has_alpha) {
    const int alpha_off = 2 * p + n_v;
    double s = 0.0;
    for (int k = 0; k < p; ++k) {
      const double a = theta[alpha_off + k];
      s += a * a;
    }
    pen += pp.alpha_lambda * s;
  }
  return pen;
}

// Accumulate the gradient of the penalty into g[0..n_theta). Only fills the
// blocks the penalty actually touches; other entries of g are left intact
// (caller must initialise them).
inline void prior_penalty_grad_add(const double* theta, int p, int n_v,
                                   bool has_alpha,
                                   const PriorParams& pp, double* g) {
  if (pp.mu_lambda > 0.0 && pp.mu_center.size() == p &&
      pp.log_sigma_center.size() == p) {
    for (int k = 0; k < p; ++k) {
      const double scale = std::exp(pp.log_sigma_center(k));
      const double d = (theta[k] - pp.mu_center(k)) / scale;
      g[k] += 2.0 * pp.mu_lambda * d / scale;
    }
  }
  if (pp.log_sigma_lambda > 0.0 && pp.log_sigma_center.size() == p) {
    for (int k = 0; k < p; ++k) {
      const double d = theta[p + k] - pp.log_sigma_center(k);
      g[p + k] += 2.0 * pp.log_sigma_lambda * d;
    }
  }
  if (pp.alpha_lambda > 0.0 && has_alpha) {
    const int alpha_off = 2 * p + n_v;
    for (int k = 0; k < p; ++k) {
      g[alpha_off + k] += 2.0 * pp.alpha_lambda * theta[alpha_off + k];
    }
  }
}

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
double loglik_niche_math_ip_weighted_eigen(
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
double loglik_niche_math_ip_weighted_grad_eigen(
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
// EXACTLY, plus the optional PriorParams penalty (ridge on log_sigma,
// optionally also on mu and alpha) that stabilises the optimization
// against the Patil & Ord (1976) sigma -> infinity drift in the pure-ML
// weighted-distribution model and against mu drifting outside the data
// cloud:
//
//   -log L = 0.5 sum_i q1(x_i)
//          - sum_i log w(x_i)
//          + n_occ * log( sum_j w(y_j) * exp(-q2(y_j)/2) )
//          + prior_penalty_value(theta, p, n_v, /*has_alpha*/false, pp)
//
// All entries of `pp` default to zero (no penalty), in which case the
// kernel reduces exactly to the unpenalised paper formula.
double loglik_niche_math_weighted_eigen(
    const double* theta, int n_theta,
    const Eigen::MatrixXd& env_occ,
    const Eigen::MatrixXd& M_den,
    const Eigen::VectorXd& w_occ,
    const Eigen::VectorXd& w_den,
    double eta,
    const PriorParams& pp);

// Hybrid analytic gradient for the penalized weighted kernel; same hybrid
// strategy as loglik_niche_math_ip_weighted_grad_eigen plus the closed-form
// gradient of every active term in PriorParams (mu / log_sigma / alpha;
// alpha is irrelevant for the Gaussian weighted family).
double loglik_niche_math_weighted_grad_eigen(
    const double* theta, int n_theta,
    const Eigen::MatrixXd& env_occ,
    const Eigen::MatrixXd& M_den,
    const Eigen::VectorXd& w_occ,
    const Eigen::VectorXd& w_den,
    double eta,
    const PriorParams& pp,
    double gradstep_rel, double gradstep_abs,
    double* g_out);

// Math-scale skew-normal presence-only negative log-likelihood. Theta layout:
//   [mu(0..p-1), log_sigma(0..p-1), v(0..p*(p-1)/2 - 1), alpha(0..p-1)]
// (the additional alpha block extends the Gaussian theta by p parameters).
//
// Density: f_{SN}(x; mu, Sigma, alpha) = 2 * phi_k(x-mu; Sigma) * Phi(z(x))
// where z(x) = sum_k alpha_k * (x_k - mu_k) / sigma_k.
// Reference: Azzalini & Capitanio (1999), J. R. Stat. Soc. Ser. B 61(3).
//
// Optionally adds a ridge penalty on alpha (well-known fix for the
// unbounded-MLE pathology of the SN direct parameterization, Azzalini 1985,
// Pewsey 2000) and/or on mu / log_sigma. mu_lambda and log_sigma_lambda are
// not commonly used for the presence-only kernel, but the same PriorParams
// is accepted for signature uniformity with the weighted variants.
double loglik_niche_math_skew_normal_eigen(
    const double* theta, int n_theta,
    const Eigen::MatrixXd& env_occ, double eta,
    const PriorParams& pp);

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
    const PriorParams& pp);

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
    const Eigen::MatrixXd& env_occ, double eta,
    const PriorParams& pp);

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
    const PriorParams& pp);

// Math-scale NCST (Non-Central Skew t) presence-only negative log-likelihood.
// Hasan & Chen (2025, arXiv:2507.10465v1): T = X / sqrt(Y/r),
// X ~ SN_k(xi, Omega, alpha), Y ~ chi^2_r. Location xi enters BEFORE
// chi-squared scaling (non-central). Same theta layout as skew-t:
//   [xi(0..p-1), log_sigma(0..p-1), v(0..p*(p-1)/2 - 1),
//    alpha(0..p-1), log_r(0)]
double loglik_niche_math_ncst_eigen(
    const double* theta, int n_theta,
    const Eigen::MatrixXd& env_occ, double eta,
    const PriorParams& pp);

// Math-scale NCST weighted (NCST density + paper Eq. 5 + ridge prior).
double loglik_niche_math_ncst_weighted_eigen(
    const double* theta, int n_theta,
    const Eigen::MatrixXd& env_occ,
    const Eigen::MatrixXd& M_den,
    const Eigen::VectorXd& w_occ,
    const Eigen::VectorXd& w_den,
    double eta,
    const PriorParams& pp);

}

#endif
