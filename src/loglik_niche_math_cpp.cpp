// src/loglik_niche_math_cpp.cpp
//
// Math-scale (theta) negative log-likelihood kernels in pure C++/Eigen.
//
// These are the kernels driven by `ucminfcpp::ucminf_xptr` via the closure
// built in `create_niche_obj_ptr` (see src/niche_obj.cpp). Compared with
// `loglik_niche_chol_cpp` and `loglik_niche_kde_bias_corrected_cpp`, they:
//
//   * accept the math-scale parameter vector `theta` directly as a
//     `const double*` (no Rcpp::NumericVector / Rcpp::NumericMatrix
//     allocations on every objective evaluation),
//   * build the C-vine correlation Cholesky factor through the Eigen-native
//     `nicher::cvine_cholesky_eigen` (no R object round-trip), and
//   * for the weighted model, take the KDE weights `w_occ` / `w_den` as
//     PRECOMPUTED inputs (KDE depends only on the environmental data, not
//     on theta, so it must NOT be recomputed each step).
//
// The hybrid analytic gradient (Option A in the implementation plan) computes
//   * grad mu, grad log_sigma  -> closed form,
//   * grad v_k                 -> central finite difference, perturbing only
//                                  the C-vine partial-correlation block.
//
// Theta layout (length n_theta = 2*p + p*(p-1)/2):
//   [mu(0..p-1),  log_sigma(0..p-1),  v(0..p*(p-1)/2 - 1)]

#include "nicher_types.h"
#include <cmath>
#include <limits>
#include <vector>

namespace nicher {

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

// Build L_corr (correlation Cholesky) and L_cov = diag(sigma) * L_corr.
// `sigma` and `v` are views into `theta`; `L_corr` and `L_cov` are reused
// scratch buffers (no per-call allocations).
static inline void build_L_cov(const double* theta, int p,
                               int n_v, double eta,
                               Eigen::VectorXd& mu_out,
                               Eigen::VectorXd& sigma_out,
                               Eigen::MatrixXd& L_corr,
                               Eigen::MatrixXd& L_cov) {
  mu_out.resize(p);
  sigma_out.resize(p);
  for (int i = 0; i < p; ++i) {
    mu_out(i) = theta[i];
    sigma_out(i) = std::exp(theta[p + i]);
  }
  Eigen::Map<const Eigen::VectorXd> v_map(theta + 2 * p, n_v);
  cvine_cholesky_eigen(v_map, p, eta, L_corr);

  L_cov = L_corr;
  for (int j = 0; j < p; ++j) {
    L_cov.row(j) *= sigma_out(j);
  }
}

// ---------------------------------------------------------------------------
// Presence-only (math scale)
// ---------------------------------------------------------------------------

double loglik_niche_math_presence_only_eigen(
    const double* theta, int n_theta,
    const Eigen::MatrixXd& env_occ, double eta) {
  const int p = env_occ.cols();
  const int n_v = p * (p - 1) / 2;
  const int expected = 2 * p + n_v;
  if (n_theta != expected) {
    Rcpp::stop("theta length mismatch (got %d, expected %d for p=%d).",
               n_theta, expected, p);
  }

  Eigen::VectorXd mu, sigma;
  Eigen::MatrixXd L_corr, L_cov;
  build_L_cov(theta, p, n_v, eta, mu, sigma, L_corr, L_cov);

  // sum_i ||L_cov^{-1} (x_i - mu)||^2
  const int n_occ = env_occ.rows();
  Eigen::MatrixXd diff(p, n_occ);
  for (int i = 0; i < n_occ; ++i) diff.col(i) = env_occ.row(i).transpose() - mu;
  Eigen::MatrixXd y = L_cov.triangularView<Eigen::Lower>().solve(diff);
  const double sum_q = y.colwise().squaredNorm().sum();

  double log_det = 0.0;
  for (int i = 0; i < p; ++i) log_det += std::log(L_cov(i, i));
  log_det *= 2.0;

  const double n = static_cast<double>(n_occ);
  double neg_log = 0.5 * n * log_det + 0.5 * sum_q;
  if (!std::isfinite(neg_log)) neg_log = OPTIM_PENALTY;
  return neg_log;
}

// ---------------------------------------------------------------------------
// Weighted (math scale, with precomputed KDE weights)
// ---------------------------------------------------------------------------

double loglik_niche_math_kde_bias_corrected_eigen(
    const double* theta, int n_theta,
    const Eigen::MatrixXd& env_occ,
    const Eigen::MatrixXd& M_den,
    const Eigen::VectorXd& w_occ,
    const Eigen::VectorXd& w_den,
    double eta) {
  const int p = env_occ.cols();
  const int n_v = p * (p - 1) / 2;
  const int expected = 2 * p + n_v;
  if (n_theta != expected) {
    Rcpp::stop("theta length mismatch (got %d, expected %d for p=%d).",
               n_theta, expected, p);
  }
  if (M_den.cols() != p) Rcpp::stop("M_den must have p columns");
  if (w_occ.size() != env_occ.rows())
    Rcpp::stop("w_occ length must equal nrow(env_occ).");
  if (w_den.size() != M_den.rows())
    Rcpp::stop("w_den length must equal nrow(M_den).");

  Eigen::VectorXd mu, sigma;
  Eigen::MatrixXd L_corr, L_cov;
  build_L_cov(theta, p, n_v, eta, mu, sigma, L_corr, L_cov);

  const int n_occ = env_occ.rows();
  const int n_den = M_den.rows();

  // Mahalanobis contributions
  Eigen::MatrixXd diff_occ(p, n_occ);
  for (int i = 0; i < n_occ; ++i) diff_occ.col(i) = env_occ.row(i).transpose() - mu;
  Eigen::MatrixXd y_occ = L_cov.triangularView<Eigen::Lower>().solve(diff_occ);
  const double sum_q1 = y_occ.colwise().squaredNorm().sum();

  Eigen::MatrixXd diff_den(p, n_den);
  for (int j = 0; j < n_den; ++j) diff_den.col(j) = M_den.row(j).transpose() - mu;
  Eigen::MatrixXd y_den = L_cov.triangularView<Eigen::Lower>().solve(diff_den);
  Eigen::ArrayXd q2 = y_den.colwise().squaredNorm().array();

  // Log-likelihood (negative). KDE weights are clamped before log() to
  // mirror the existing C++ kernel (loglik_niche_kde_bias_corrected_cpp).
  Eigen::ArrayXd log_w_occ = w_occ.array().max(MIN_KDE_WEIGHT).log();
  Eigen::ArrayXd log_w_den = w_den.array().max(MIN_KDE_WEIGHT).log();
  Eigen::ArrayXd a = -0.5 * q2 - log_w_den;

  const double max_a = a.maxCoeff();
  const double sum_exp = (a - max_a).exp().sum();
  const double log_sum_exp = max_a + std::log(sum_exp);

  double neg_log = 0.5 * sum_q1
                 + log_w_occ.sum()
                 + static_cast<double>(n_occ) * log_sum_exp;
  if (!std::isfinite(neg_log)) neg_log = OPTIM_PENALTY;
  return neg_log;
}

// ---------------------------------------------------------------------------
// Hybrid analytic gradient (Option A)
//
// f(theta) = 0.5 sum_i ||y_i||^2  +  sum_i log w(x_i)  +  n_occ * Z
//   y_i = L^{-1} (x_i - mu),     z_j = L^{-1} (m_j - mu)
//   Z   = log sum_j exp(a_j),    a_j = -0.5 ||z_j||^2 - log w(m_j)
//   pi_j = softmax(a_j), summing to 1.
//
// Because the KDE weights w(.) depend only on env data (not on theta),
// d log w / d theta = 0 -> the gradient reduces to the standard Gaussian
// gradient with the M-side reweighted by pi.
//
// With L = diag(sigma) * L_corr, change variables u = (.) ./ sigma. Then
//   y = L_corr^{-1} u,   v := L_corr^{-T} y,   v = M_corr u  with
//   M_corr = (L_corr L_corr^T)^{-1}.
//
//   d/d s_k (0.5 ||y||^2) = - u_k * v_k    (because d u / d s_k = -u_k e_k)
//
// So with U_occ = diag(1/sigma) (X_occ - mu)^T (p x n_occ),
//        Y_occ = L_corr^{-1} U_occ, V_occ = L_corr^{-T} Y_occ:
//   d f_pres / d s_k = - sum_i (U_occ ⊙ V_occ)_{k, i}
// and analogously for the M-side, weighted by pi_j.
//
// d f / d mu = Sigma^{-1} [ -sum_i (x_i - mu) + n_occ sum_j pi_j (m_j - mu) ]
//            = L^{-T} L^{-1} [ -sum_i (x_i - mu) + n_occ * (M_den^T pi - mu) ]
//
// The v-block uses central FD by perturbing theta[2p + k] in place; the M-side
// KDE is fixed so each sub-evaluation is one full kernel call.
// ---------------------------------------------------------------------------

double loglik_niche_math_kde_bias_corrected_grad_eigen(
    const double* theta, int n_theta,
    const Eigen::MatrixXd& env_occ,
    const Eigen::MatrixXd& M_den,
    const Eigen::VectorXd& w_occ,
    const Eigen::VectorXd& w_den,
    double eta,
    double gradstep_rel, double gradstep_abs,
    double* g_out) {
  const int p = env_occ.cols();
  const int n_v = p * (p - 1) / 2;
  const int expected = 2 * p + n_v;
  if (n_theta != expected) {
    Rcpp::stop("theta length mismatch (got %d, expected %d for p=%d).",
               n_theta, expected, p);
  }

  // --- 1. Build mu, sigma, L_corr, L_cov (shared with f-only path) -------
  Eigen::VectorXd mu, sigma;
  Eigen::MatrixXd L_corr, L_cov;
  build_L_cov(theta, p, n_v, eta, mu, sigma, L_corr, L_cov);

  const int n_occ = env_occ.rows();
  const int n_den = M_den.rows();

  // --- 2. Forward solves: U, Y, V on each side ---------------------------
  // U_occ(k, i) = (env_occ(i, k) - mu_k) / sigma_k
  Eigen::MatrixXd U_occ(p, n_occ);
  for (int i = 0; i < n_occ; ++i) {
    U_occ.col(i) = (env_occ.row(i).transpose() - mu).cwiseQuotient(sigma);
  }
  Eigen::MatrixXd U_den(p, n_den);
  for (int j = 0; j < n_den; ++j) {
    U_den.col(j) = (M_den.row(j).transpose() - mu).cwiseQuotient(sigma);
  }

  Eigen::MatrixXd Y_occ = L_corr.triangularView<Eigen::Lower>().solve(U_occ);
  Eigen::MatrixXd Y_den = L_corr.triangularView<Eigen::Lower>().solve(U_den);

  Eigen::MatrixXd V_occ = L_corr.transpose().triangularView<Eigen::Upper>().solve(Y_occ);
  Eigen::MatrixXd V_den = L_corr.transpose().triangularView<Eigen::Upper>().solve(Y_den);

  // --- 3. f value (matches loglik_niche_math_kde_bias_corrected_eigen exactly) -----
  const double sum_q1 = Y_occ.colwise().squaredNorm().sum();
  Eigen::ArrayXd q2 = Y_den.colwise().squaredNorm().array();

  Eigen::ArrayXd log_w_occ = w_occ.array().max(MIN_KDE_WEIGHT).log();
  Eigen::ArrayXd log_w_den = w_den.array().max(MIN_KDE_WEIGHT).log();
  Eigen::ArrayXd a = -0.5 * q2 - log_w_den;

  const double max_a = a.maxCoeff();
  const Eigen::ArrayXd ea = (a - max_a).exp();
  const double sum_exp = ea.sum();
  const double log_sum_exp = max_a + std::log(sum_exp);

  double f = 0.5 * sum_q1
           + log_w_occ.sum()
           + static_cast<double>(n_occ) * log_sum_exp;
  if (!std::isfinite(f)) {
    f = OPTIM_PENALTY;
    std::fill(g_out, g_out + n_theta, 0.0);
    return f;
  }

  // Softmax weights (sum to 1)
  Eigen::VectorXd pi = (ea / sum_exp).matrix();

  // --- 4. Analytic grad w.r.t. mu ----------------------------------------
  //   sum_d_occ = sum_i (x_i - mu)   (p)
  //   weighted_e = M_den^T pi - mu   (p)
  Eigen::VectorXd sum_d_occ = env_occ.colwise().sum().transpose()
                               - static_cast<double>(n_occ) * mu;
  Eigen::VectorXd weighted_e = M_den.transpose() * pi - mu;
  Eigen::VectorXd rhs = -sum_d_occ + static_cast<double>(n_occ) * weighted_e;

  Eigen::VectorXd grad_mu = L_cov.triangularView<Eigen::Lower>().solve(rhs);
  grad_mu = L_cov.transpose().triangularView<Eigen::Upper>().solve(grad_mu);

  // --- 5. Analytic grad w.r.t. log_sigma ---------------------------------
  //   d f / d s_k = -sum_i U_occ(k,i) V_occ(k,i)
  //                + n_occ sum_j pi_j U_den(k,j) V_den(k,j)
  Eigen::VectorXd grad_s(p);
  Eigen::ArrayXXd UV_occ = U_occ.array() * V_occ.array();        // p x n_occ
  Eigen::ArrayXXd UV_den = U_den.array() * V_den.array();        // p x n_den
  Eigen::VectorXd term_pres = UV_occ.matrix().rowwise().sum();   // p
  Eigen::VectorXd term_den  = UV_den.matrix() * pi;              // p
  grad_s = -term_pres + static_cast<double>(n_occ) * term_den;

  // --- 6. FD grad w.r.t. v block (central differences) -------------------
  std::vector<double> theta_pert(theta, theta + n_theta);
  Eigen::VectorXd grad_v(n_v);
  for (int k = 0; k < n_v; ++k) {
    const int idx = 2 * p + k;
    const double xk = theta_pert[idx];
    const double dx = std::abs(xk) * gradstep_rel + gradstep_abs;

    theta_pert[idx] = xk + dx;
    const double f_plus = loglik_niche_math_kde_bias_corrected_eigen(
        theta_pert.data(), n_theta, env_occ, M_den, w_occ, w_den, eta);

    theta_pert[idx] = xk - dx;
    const double f_minus = loglik_niche_math_kde_bias_corrected_eigen(
        theta_pert.data(), n_theta, env_occ, M_den, w_occ, w_den, eta);

    theta_pert[idx] = xk;
    grad_v(k) = (f_plus - f_minus) / (2.0 * dx);
  }

  // --- 7. Splice grad blocks into g_out ----------------------------------
  for (int i = 0; i < p; ++i) g_out[i] = grad_mu(i);
  for (int i = 0; i < p; ++i) g_out[p + i] = grad_s(i);
  for (int k = 0; k < n_v; ++k) g_out[2 * p + k] = grad_v(k);

  // Sanitize against NaN/Inf in any grad component
  for (int i = 0; i < n_theta; ++i) {
    if (!std::isfinite(g_out[i])) g_out[i] = 0.0;
  }
  return f;
}

// ---------------------------------------------------------------------------
// Penalized weighted (math scale): paper Eq. 5 + ridge prior on log_sigma
//
// f(theta) = 0.5 sum_i ||y_i||^2  -  sum_i log w(x_i)  +  n_occ * Z
//          + lambda * sum_k (log_sigma_k - log_sigma_center_k)^2
//
//   y_i = L^{-1} (x_i - mu),     z_j = L^{-1} (m_j - mu)
//   Z   = log sum_j exp(a_j),    a_j = log w(m_j) - 0.5 ||z_j||^2
//   pi_j = softmax(a_j), summing to 1.
//
// This is the formula written in Jimenez & Soberon 2022 (Ecological
// Modelling 438:109982), Eq. 5, with a weakly-informative ridge penalty
// added on the log-standard-deviations to prevent the well-known
// Patil & Ord (1976) sigma -> infinity drift in the pure-ML solution.
// The prior is centered at log_sigma_center (typically log(sd(env_occ))
// computed by the R caller); lambda controls the strength.
// ---------------------------------------------------------------------------

double loglik_niche_math_weighted_eigen(
    const double* theta, int n_theta,
    const Eigen::MatrixXd& env_occ,
    const Eigen::MatrixXd& M_den,
    const Eigen::VectorXd& w_occ,
    const Eigen::VectorXd& w_den,
    double eta,
    const Eigen::VectorXd& prior_log_sigma_center,
    double prior_log_sigma_lambda) {
  const int p = env_occ.cols();
  const int n_v = p * (p - 1) / 2;
  const int expected = 2 * p + n_v;
  if (n_theta != expected) {
    Rcpp::stop("theta length mismatch (got %d, expected %d for p=%d).",
               n_theta, expected, p);
  }
  if (M_den.cols() != p) Rcpp::stop("M_den must have p columns");
  if (w_occ.size() != env_occ.rows())
    Rcpp::stop("w_occ length must equal nrow(env_occ).");
  if (w_den.size() != M_den.rows())
    Rcpp::stop("w_den length must equal nrow(M_den).");
  if (prior_log_sigma_center.size() != p)
    Rcpp::stop("prior_log_sigma_center length must equal p (got %d, expected %d).",
               (int)prior_log_sigma_center.size(), p);
  if (!(prior_log_sigma_lambda >= 0.0))
    Rcpp::stop("prior_log_sigma_lambda must be non-negative.");

  Eigen::VectorXd mu, sigma;
  Eigen::MatrixXd L_corr, L_cov;
  build_L_cov(theta, p, n_v, eta, mu, sigma, L_corr, L_cov);

  const int n_occ = env_occ.rows();
  const int n_den = M_den.rows();

  Eigen::MatrixXd diff_occ(p, n_occ);
  for (int i = 0; i < n_occ; ++i) diff_occ.col(i) = env_occ.row(i).transpose() - mu;
  Eigen::MatrixXd y_occ = L_cov.triangularView<Eigen::Lower>().solve(diff_occ);
  const double sum_q1 = y_occ.colwise().squaredNorm().sum();

  Eigen::MatrixXd diff_den(p, n_den);
  for (int j = 0; j < n_den; ++j) diff_den.col(j) = M_den.row(j).transpose() - mu;
  Eigen::MatrixXd y_den = L_cov.triangularView<Eigen::Lower>().solve(diff_den);
  Eigen::ArrayXd q2 = y_den.colwise().squaredNorm().array();

  // Paper Eq. 5: a_j = log w(y_j) - 0.5 q2_j, and the outer sum has
  //   - sum_i log w(x_i) (NOT + as in the legacy `weighted` formula).
  Eigen::ArrayXd log_w_occ = w_occ.array().max(MIN_KDE_WEIGHT).log();
  Eigen::ArrayXd log_w_den = w_den.array().max(MIN_KDE_WEIGHT).log();
  Eigen::ArrayXd a = log_w_den - 0.5 * q2;

  const double max_a = a.maxCoeff();
  const double sum_exp = (a - max_a).exp().sum();
  const double log_sum_exp = max_a + std::log(sum_exp);

  // Ridge penalty on log_sigma. theta layout has log_sigma at offsets [p..2p-1].
  double ridge = 0.0;
  if (prior_log_sigma_lambda > 0.0) {
    for (int k = 0; k < p; ++k) {
      const double d = theta[p + k] - prior_log_sigma_center(k);
      ridge += d * d;
    }
    ridge *= prior_log_sigma_lambda;
  }

  double neg_log = 0.5 * sum_q1
                 - log_w_occ.sum()
                 + static_cast<double>(n_occ) * log_sum_exp
                 + ridge;
  if (!std::isfinite(neg_log)) neg_log = OPTIM_PENALTY;
  return neg_log;
}

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
    double* g_out) {
  const int p = env_occ.cols();
  const int n_v = p * (p - 1) / 2;
  const int expected = 2 * p + n_v;
  if (n_theta != expected) {
    Rcpp::stop("theta length mismatch (got %d, expected %d for p=%d).",
               n_theta, expected, p);
  }
  if (prior_log_sigma_center.size() != p)
    Rcpp::stop("prior_log_sigma_center length must equal p (got %d, expected %d).",
               (int)prior_log_sigma_center.size(), p);
  if (!(prior_log_sigma_lambda >= 0.0))
    Rcpp::stop("prior_log_sigma_lambda must be non-negative.");

  Eigen::VectorXd mu, sigma;
  Eigen::MatrixXd L_corr, L_cov;
  build_L_cov(theta, p, n_v, eta, mu, sigma, L_corr, L_cov);

  const int n_occ = env_occ.rows();
  const int n_den = M_den.rows();

  Eigen::MatrixXd U_occ(p, n_occ);
  for (int i = 0; i < n_occ; ++i)
    U_occ.col(i) = (env_occ.row(i).transpose() - mu).cwiseQuotient(sigma);
  Eigen::MatrixXd U_den(p, n_den);
  for (int j = 0; j < n_den; ++j)
    U_den.col(j) = (M_den.row(j).transpose() - mu).cwiseQuotient(sigma);

  Eigen::MatrixXd Y_occ = L_corr.triangularView<Eigen::Lower>().solve(U_occ);
  Eigen::MatrixXd Y_den = L_corr.triangularView<Eigen::Lower>().solve(U_den);

  Eigen::MatrixXd V_occ = L_corr.transpose().triangularView<Eigen::Upper>().solve(Y_occ);
  Eigen::MatrixXd V_den = L_corr.transpose().triangularView<Eigen::Upper>().solve(Y_den);

  const double sum_q1 = Y_occ.colwise().squaredNorm().sum();
  Eigen::ArrayXd q2 = Y_den.colwise().squaredNorm().array();

  Eigen::ArrayXd log_w_occ = w_occ.array().max(MIN_KDE_WEIGHT).log();
  Eigen::ArrayXd log_w_den = w_den.array().max(MIN_KDE_WEIGHT).log();
  Eigen::ArrayXd a = log_w_den - 0.5 * q2;

  const double max_a = a.maxCoeff();
  const Eigen::ArrayXd ea = (a - max_a).exp();
  const double sum_exp = ea.sum();
  const double log_sum_exp = max_a + std::log(sum_exp);

  double ridge = 0.0;
  if (prior_log_sigma_lambda > 0.0) {
    for (int k = 0; k < p; ++k) {
      const double d = theta[p + k] - prior_log_sigma_center(k);
      ridge += d * d;
    }
    ridge *= prior_log_sigma_lambda;
  }

  double f = 0.5 * sum_q1
           - log_w_occ.sum()
           + static_cast<double>(n_occ) * log_sum_exp
           + ridge;
  if (!std::isfinite(f)) {
    f = OPTIM_PENALTY;
    std::fill(g_out, g_out + n_theta, 0.0);
    return f;
  }

  // Softmax weights pi_j (sum to 1)
  Eigen::VectorXd pi = (ea / sum_exp).matrix();

  // Gradient w.r.t. mu (log w terms are constant in theta, so identical
  // structure to loglik_niche_math_kde_bias_corrected_grad_eigen).
  Eigen::VectorXd sum_d_occ = env_occ.colwise().sum().transpose()
                              - static_cast<double>(n_occ) * mu;
  Eigen::VectorXd weighted_e = M_den.transpose() * pi - mu;
  Eigen::VectorXd rhs = -sum_d_occ + static_cast<double>(n_occ) * weighted_e;

  Eigen::VectorXd grad_mu = L_cov.triangularView<Eigen::Lower>().solve(rhs);
  grad_mu = L_cov.transpose().triangularView<Eigen::Upper>().solve(grad_mu);

  // Gradient w.r.t. log_sigma (Gaussian part, then add ridge derivative).
  Eigen::ArrayXXd UV_occ = U_occ.array() * V_occ.array();
  Eigen::ArrayXXd UV_den = U_den.array() * V_den.array();
  Eigen::VectorXd term_pres = UV_occ.matrix().rowwise().sum();
  Eigen::VectorXd term_den  = UV_den.matrix() * pi;
  Eigen::VectorXd grad_s = -term_pres + static_cast<double>(n_occ) * term_den;
  if (prior_log_sigma_lambda > 0.0) {
    for (int k = 0; k < p; ++k) {
      grad_s(k) += 2.0 * prior_log_sigma_lambda *
                   (theta[p + k] - prior_log_sigma_center(k));
    }
  }

  // Gradient w.r.t. v block via central FD (same recipe as the unpenalized
  // weighted gradient; the ridge does not depend on v).
  std::vector<double> theta_pert(theta, theta + n_theta);
  Eigen::VectorXd grad_v(n_v);
  for (int k = 0; k < n_v; ++k) {
    const int idx = 2 * p + k;
    const double xk = theta_pert[idx];
    const double dx = std::abs(xk) * gradstep_rel + gradstep_abs;

    theta_pert[idx] = xk + dx;
    const double f_plus = loglik_niche_math_weighted_eigen(
        theta_pert.data(), n_theta, env_occ, M_den, w_occ, w_den, eta,
        prior_log_sigma_center, prior_log_sigma_lambda);

    theta_pert[idx] = xk - dx;
    const double f_minus = loglik_niche_math_weighted_eigen(
        theta_pert.data(), n_theta, env_occ, M_den, w_occ, w_den, eta,
        prior_log_sigma_center, prior_log_sigma_lambda);

    theta_pert[idx] = xk;
    grad_v(k) = (f_plus - f_minus) / (2.0 * dx);
  }

  for (int i = 0; i < p; ++i) g_out[i] = grad_mu(i);
  for (int i = 0; i < p; ++i) g_out[p + i] = grad_s(i);
  for (int k = 0; k < n_v; ++k) g_out[2 * p + k] = grad_v(k);

  for (int i = 0; i < n_theta; ++i) {
    if (!std::isfinite(g_out[i])) g_out[i] = 0.0;
  }
  return f;
}

// ---------------------------------------------------------------------------
// Skew-normal helpers
// ---------------------------------------------------------------------------
//
// Theta layout for skew_normal (length n_theta = 3*p + p*(p-1)/2):
//   [mu(0..p-1),  log_sigma(0..p-1),  v(0..p*(p-1)/2 - 1),  alpha(0..p-1)]
//
// Density: f_{SN}(x; mu, Sigma, alpha) = 2 * phi_k(x-mu; Sigma) * Phi(z(x))
// where Phi is the univariate standard-normal CDF and
//   z(x) = sum_k alpha_k * (x_k - mu_k) / sigma_k.
// (Azzalini & Capitanio 1999, J. R. Stat. Soc. Ser. B 61(3): 579-602.)

// log Phi(z), numerically stable (uses R's pnorm with log_p=TRUE)
static inline double log_pnorm(double z) {
  return R::pnorm(z, 0.0, 1.0, /*lower_tail=*/1, /*log_p=*/1);
}

// Compute z_i = sum_k alpha_k * (x_ki - mu_k) / sigma_k for every column of
// (X - mu), where X has p rows and N columns.
static inline Eigen::VectorXd skew_z(const Eigen::MatrixXd& diff,
                                     const Eigen::VectorXd& sigma,
                                     const Eigen::VectorXd& alpha) {
  // alpha_over_sigma_k = alpha_k / sigma_k
  Eigen::VectorXd a_over_s = alpha.array() / sigma.array();
  // z_i = a_over_s^T diff_i
  return diff.transpose() * a_over_s;
}

// Sum of log Phi(z_i) over the columns of `diff`.
static inline double sum_log_pnorm(const Eigen::MatrixXd& diff,
                                   const Eigen::VectorXd& sigma,
                                   const Eigen::VectorXd& alpha) {
  Eigen::VectorXd z = skew_z(diff, sigma, alpha);
  double s = 0.0;
  for (int i = 0; i < z.size(); ++i) s += log_pnorm(z(i));
  return s;
}

// ---------------------------------------------------------------------------
// Skew-normal presence-only (math scale)
// ---------------------------------------------------------------------------

double loglik_niche_math_skew_normal_eigen(
    const double* theta, int n_theta,
    const Eigen::MatrixXd& env_occ, double eta) {
  const int p = env_occ.cols();
  const int n_v = p * (p - 1) / 2;
  const int expected = 3 * p + n_v;
  if (n_theta != expected) {
    Rcpp::stop("theta length mismatch (got %d, expected %d for p=%d skew-normal).",
               n_theta, expected, p);
  }

  Eigen::VectorXd mu, sigma;
  Eigen::MatrixXd L_corr, L_cov;
  build_L_cov(theta, p, n_v, eta, mu, sigma, L_corr, L_cov);

  Eigen::Map<const Eigen::VectorXd> alpha(theta + 2 * p + n_v, p);

  const int n_occ = env_occ.rows();
  Eigen::MatrixXd diff(p, n_occ);
  for (int i = 0; i < n_occ; ++i) diff.col(i) = env_occ.row(i).transpose() - mu;

  Eigen::MatrixXd y = L_cov.triangularView<Eigen::Lower>().solve(diff);
  const double sum_q = y.colwise().squaredNorm().sum();

  double log_det = 0.0;
  for (int i = 0; i < p; ++i) log_det += std::log(L_cov(i, i));
  log_det *= 2.0;

  const double sum_log_phi = sum_log_pnorm(diff, sigma, alpha);

  // Constants n*log(2) and -n*k/2*log(2*pi) drop out (no influence on the
  // optimum and they would not match the Gaussian PO kernel which also
  // drops the (2*pi) constant).
  const double n = static_cast<double>(n_occ);
  double neg_log = 0.5 * n * log_det + 0.5 * sum_q - sum_log_phi;
  if (!std::isfinite(neg_log)) neg_log = OPTIM_PENALTY;
  return neg_log;
}

// ---------------------------------------------------------------------------
// Skew-normal weighted (math scale, with precomputed KDE weights, paper Eq. 5
// + ridge prior on log_sigma)
// ---------------------------------------------------------------------------

double loglik_niche_math_skew_normal_weighted_eigen(
    const double* theta, int n_theta,
    const Eigen::MatrixXd& env_occ,
    const Eigen::MatrixXd& M_den,
    const Eigen::VectorXd& w_occ,
    const Eigen::VectorXd& w_den,
    double eta,
    const Eigen::VectorXd& prior_log_sigma_center,
    double prior_log_sigma_lambda) {
  const int p = env_occ.cols();
  const int n_v = p * (p - 1) / 2;
  const int expected = 3 * p + n_v;
  if (n_theta != expected) {
    Rcpp::stop("theta length mismatch (got %d, expected %d for p=%d skew-normal).",
               n_theta, expected, p);
  }
  if (M_den.cols() != p) Rcpp::stop("M_den must have p columns");
  if (w_occ.size() != env_occ.rows())
    Rcpp::stop("w_occ length must equal nrow(env_occ).");
  if (w_den.size() != M_den.rows())
    Rcpp::stop("w_den length must equal nrow(M_den).");
  if (prior_log_sigma_center.size() != p)
    Rcpp::stop("prior_log_sigma_center length must equal p (got %d, expected %d).",
               (int)prior_log_sigma_center.size(), p);
  if (!(prior_log_sigma_lambda >= 0.0))
    Rcpp::stop("prior_log_sigma_lambda must be non-negative.");

  Eigen::VectorXd mu, sigma;
  Eigen::MatrixXd L_corr, L_cov;
  build_L_cov(theta, p, n_v, eta, mu, sigma, L_corr, L_cov);

  Eigen::Map<const Eigen::VectorXd> alpha(theta + 2 * p + n_v, p);

  const int n_occ = env_occ.rows();
  const int n_den = M_den.rows();

  Eigen::MatrixXd diff_occ(p, n_occ);
  for (int i = 0; i < n_occ; ++i) diff_occ.col(i) = env_occ.row(i).transpose() - mu;
  Eigen::MatrixXd y_occ = L_cov.triangularView<Eigen::Lower>().solve(diff_occ);
  const double sum_q1 = y_occ.colwise().squaredNorm().sum();
  const double sum_log_phi_occ = sum_log_pnorm(diff_occ, sigma, alpha);

  Eigen::MatrixXd diff_den(p, n_den);
  for (int j = 0; j < n_den; ++j) diff_den.col(j) = M_den.row(j).transpose() - mu;
  Eigen::MatrixXd y_den = L_cov.triangularView<Eigen::Lower>().solve(diff_den);
  Eigen::ArrayXd q2 = y_den.colwise().squaredNorm().array();

  // log Phi(z) for denominator points
  Eigen::VectorXd z_den = skew_z(diff_den, sigma, alpha);
  Eigen::ArrayXd log_phi_den(n_den);
  for (int j = 0; j < n_den; ++j) log_phi_den(j) = log_pnorm(z_den(j));

  Eigen::ArrayXd log_w_occ = w_occ.array().max(MIN_KDE_WEIGHT).log();
  Eigen::ArrayXd log_w_den = w_den.array().max(MIN_KDE_WEIGHT).log();

  // Inner integrand (log scale): log w_den + log Phi(z_den) - 0.5 q2
  Eigen::ArrayXd a = log_w_den + log_phi_den - 0.5 * q2;

  const double max_a = a.maxCoeff();
  const double sum_exp = (a - max_a).exp().sum();
  const double log_sum_exp = max_a + std::log(sum_exp);

  // Ridge penalty on log_sigma. theta layout has log_sigma at offsets [p..2p-1].
  double ridge = 0.0;
  if (prior_log_sigma_lambda > 0.0) {
    for (int k = 0; k < p; ++k) {
      const double d = theta[p + k] - prior_log_sigma_center(k);
      ridge += d * d;
    }
    ridge *= prior_log_sigma_lambda;
  }

  double neg_log = 0.5 * sum_q1
                 - sum_log_phi_occ
                 - log_w_occ.sum()
                 + static_cast<double>(n_occ) * log_sum_exp
                 + ridge;
  if (!std::isfinite(neg_log)) neg_log = OPTIM_PENALTY;
  return neg_log;
}

} // namespace nicher

// ===========================================================================
// Rcpp-exported wrappers (used for parity testing from R; the optimizer hot
// path goes through the XPtr in `niche_obj.cpp` which calls the namespaced
// helpers above directly).
// ===========================================================================

using namespace Rcpp;

// [[Rcpp::export]]
double loglik_niche_math_presence_only_cpp(
    NumericVector theta, NumericMatrix env_occ, double eta = 1.0) {
  Eigen::Map<Eigen::MatrixXd> occ_map(
      Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(env_occ));
  Eigen::MatrixXd occ = occ_map;
  return nicher::loglik_niche_math_presence_only_eigen(
      &theta[0], theta.size(), occ, eta);
}

// [[Rcpp::export]]
double loglik_niche_math_kde_bias_corrected_cpp(
    NumericVector theta,
    NumericMatrix env_occ,
    NumericMatrix M_den,
    NumericVector w_occ,
    NumericVector w_den,
    double eta = 1.0) {
  Eigen::Map<Eigen::MatrixXd> occ_map(
      Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(env_occ));
  Eigen::Map<Eigen::MatrixXd> mden_map(
      Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(M_den));
  Eigen::Map<Eigen::VectorXd> wocc_map(
      Rcpp::as<Eigen::Map<Eigen::VectorXd>>(w_occ));
  Eigen::Map<Eigen::VectorXd> wden_map(
      Rcpp::as<Eigen::Map<Eigen::VectorXd>>(w_den));

  Eigen::MatrixXd occ = occ_map, mden = mden_map;
  Eigen::VectorXd wocc = wocc_map, wden = wden_map;

  return nicher::loglik_niche_math_kde_bias_corrected_eigen(
      &theta[0], theta.size(), occ, mden, wocc, wden, eta);
}

// [[Rcpp::export]]
List loglik_niche_math_kde_bias_corrected_grad_cpp(
    NumericVector theta,
    NumericMatrix env_occ,
    NumericMatrix M_den,
    NumericVector w_occ,
    NumericVector w_den,
    double eta = 1.0,
    double gradstep_rel = 1e-6,
    double gradstep_abs = 1e-8) {
  Eigen::Map<Eigen::MatrixXd> occ_map(
      Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(env_occ));
  Eigen::Map<Eigen::MatrixXd> mden_map(
      Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(M_den));
  Eigen::Map<Eigen::VectorXd> wocc_map(
      Rcpp::as<Eigen::Map<Eigen::VectorXd>>(w_occ));
  Eigen::Map<Eigen::VectorXd> wden_map(
      Rcpp::as<Eigen::Map<Eigen::VectorXd>>(w_den));

  Eigen::MatrixXd occ = occ_map, mden = mden_map;
  Eigen::VectorXd wocc = wocc_map, wden = wden_map;

  NumericVector g(theta.size());
  double f = nicher::loglik_niche_math_kde_bias_corrected_grad_eigen(
      &theta[0], theta.size(), occ, mden, wocc, wden, eta,
      gradstep_rel, gradstep_abs, &g[0]);

  return List::create(_["value"] = f, _["gradient"] = g);
}

// [[Rcpp::export]]
double loglik_niche_math_weighted_cpp(
    NumericVector theta,
    NumericMatrix env_occ,
    NumericMatrix M_den,
    NumericVector w_occ,
    NumericVector w_den,
    NumericVector prior_log_sigma_center,
    double prior_log_sigma_lambda,
    double eta = 1.0) {
  Eigen::Map<Eigen::MatrixXd> occ_map(
      Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(env_occ));
  Eigen::Map<Eigen::MatrixXd> mden_map(
      Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(M_den));
  Eigen::Map<Eigen::VectorXd> wocc_map(
      Rcpp::as<Eigen::Map<Eigen::VectorXd>>(w_occ));
  Eigen::Map<Eigen::VectorXd> wden_map(
      Rcpp::as<Eigen::Map<Eigen::VectorXd>>(w_den));
  Eigen::Map<Eigen::VectorXd> pc_map(
      Rcpp::as<Eigen::Map<Eigen::VectorXd>>(prior_log_sigma_center));

  Eigen::MatrixXd occ = occ_map, mden = mden_map;
  Eigen::VectorXd wocc = wocc_map, wden = wden_map, pc = pc_map;

  return nicher::loglik_niche_math_weighted_eigen(
      &theta[0], theta.size(), occ, mden, wocc, wden, eta,
      pc, prior_log_sigma_lambda);
}

// [[Rcpp::export]]
List loglik_niche_math_weighted_grad_cpp(
    NumericVector theta,
    NumericMatrix env_occ,
    NumericMatrix M_den,
    NumericVector w_occ,
    NumericVector w_den,
    NumericVector prior_log_sigma_center,
    double prior_log_sigma_lambda,
    double eta = 1.0,
    double gradstep_rel = 1e-6,
    double gradstep_abs = 1e-8) {
  Eigen::Map<Eigen::MatrixXd> occ_map(
      Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(env_occ));
  Eigen::Map<Eigen::MatrixXd> mden_map(
      Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(M_den));
  Eigen::Map<Eigen::VectorXd> wocc_map(
      Rcpp::as<Eigen::Map<Eigen::VectorXd>>(w_occ));
  Eigen::Map<Eigen::VectorXd> wden_map(
      Rcpp::as<Eigen::Map<Eigen::VectorXd>>(w_den));
  Eigen::Map<Eigen::VectorXd> pc_map(
      Rcpp::as<Eigen::Map<Eigen::VectorXd>>(prior_log_sigma_center));

  Eigen::MatrixXd occ = occ_map, mden = mden_map;
  Eigen::VectorXd wocc = wocc_map, wden = wden_map, pc = pc_map;

  NumericVector g(theta.size());
  double f = nicher::loglik_niche_math_weighted_grad_eigen(
      &theta[0], theta.size(), occ, mden, wocc, wden, eta,
      pc, prior_log_sigma_lambda,
      gradstep_rel, gradstep_abs, &g[0]);

  return List::create(_["value"] = f, _["gradient"] = g);
}

// [[Rcpp::export]]
double loglik_niche_math_skew_normal_cpp(
    NumericVector theta, NumericMatrix env_occ, double eta = 1.0) {
  Eigen::Map<Eigen::MatrixXd> occ_map(
      Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(env_occ));
  Eigen::MatrixXd occ = occ_map;
  return nicher::loglik_niche_math_skew_normal_eigen(
      &theta[0], theta.size(), occ, eta);
}

// [[Rcpp::export]]
double loglik_niche_math_skew_normal_weighted_cpp(
    NumericVector theta,
    NumericMatrix env_occ,
    NumericMatrix M_den,
    NumericVector w_occ,
    NumericVector w_den,
    NumericVector prior_log_sigma_center,
    double prior_log_sigma_lambda,
    double eta = 1.0) {
  Eigen::Map<Eigen::MatrixXd> occ_map(
      Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(env_occ));
  Eigen::Map<Eigen::MatrixXd> mden_map(
      Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(M_den));
  Eigen::Map<Eigen::VectorXd> wocc_map(
      Rcpp::as<Eigen::Map<Eigen::VectorXd>>(w_occ));
  Eigen::Map<Eigen::VectorXd> wden_map(
      Rcpp::as<Eigen::Map<Eigen::VectorXd>>(w_den));
  Eigen::Map<Eigen::VectorXd> pc_map(
      Rcpp::as<Eigen::Map<Eigen::VectorXd>>(prior_log_sigma_center));

  Eigen::MatrixXd occ = occ_map, mden = mden_map;
  Eigen::VectorXd wocc = wocc_map, wden = wden_map, pc = pc_map;

  return nicher::loglik_niche_math_skew_normal_weighted_eigen(
      &theta[0], theta.size(), occ, mden, wocc, wden, eta,
      pc, prior_log_sigma_lambda);
}
