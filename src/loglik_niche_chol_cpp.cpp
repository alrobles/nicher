#include "nicher_types.h"
#include <vector>
#include <limits>

using namespace Rcpp;

// [[Rcpp::export]]
double loglik_niche_chol_cpp(NumericVector mu,
                             NumericMatrix L,          // lower Cholesky factor of Sigma
                             NumericMatrix env_occ,
                             NumericMatrix env_m) {
  int p = L.nrow();
  int n_occ = env_occ.nrow();
  int n_m = env_m.nrow();
  
  // Input validation
  if (L.ncol() != p) stop("L must be square");
  if (env_occ.ncol() != p) stop("env_occ must have p columns");
  if (env_m.ncol() != p) stop("env_m must have p columns");
  if (mu.size() != p) stop("mu must have length p");
  if (n_occ <= 0 || n_m <= 0) stop("env_occ/env_m must have n>0");
  
  // Map to Eigen (zero-copy)
  Eigen::Map<Eigen::VectorXd> mu_eig(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(mu));
  Eigen::Map<Eigen::MatrixXd> L_eig(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(L));
  Eigen::Map<Eigen::MatrixXd> occ_eig(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(env_occ));
  Eigen::Map<Eigen::MatrixXd> m_eig(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(env_m));
  
  // --- 2D specialization (maximum performance) ---
  if (p == 2) {
    double L11 = L_eig(0,0);
    double L21 = L_eig(1,0);
    double L22 = L_eig(1,1);
    double mu1 = mu_eig(0);
    double mu2 = mu_eig(1);
    double invL11 = 1.0 / L11;
    double invL22 = 1.0 / L22;
    
    double sum_q1 = 0.0;
    for (int i = 0; i < n_occ; ++i) {
      double d1 = occ_eig(i,0) - mu1;
      double d2 = occ_eig(i,1) - mu2;
      double y1 = d1 * invL11;
      double y2 = (d2 - L21 * y1) * invL22;
      sum_q1 += y1 * y1 + y2 * y2;
    }
    
    std::vector<double> a(n_m);
    double max_a = -std::numeric_limits<double>::infinity();
    for (int j = 0; j < n_m; ++j) {
      double d1 = m_eig(j,0) - mu1;
      double d2 = m_eig(j,1) - mu2;
      double y1 = d1 * invL11;
      double y2 = (d2 - L21 * y1) * invL22;
      double d2_val = y1 * y1 + y2 * y2;
      double aj = -0.5 * d2_val;
      a[j] = aj;
      if (aj > max_a) max_a = aj;
    }
    
    double sum_exp = 0.0;
    for (int j = 0; j < n_m; ++j) {
      sum_exp += std::exp(a[j] - max_a);
    }
    double log_sum_exp = max_a + std::log(sum_exp);
    double neg_log = 0.5 * sum_q1 + static_cast<double>(n_occ) * log_sum_exp;
    
    if (!std::isfinite(neg_log)) {
      neg_log = nicher::OPTIM_PENALTY;
    }
    
    return neg_log;
  }
  
  // --- Generic version for p > 2 (vectorized with Eigen) ---
  
  // Compute differences: (each column = point - mu)
  // Use arrays for broadcasting then convert to matrix to solve triangular system
  Eigen::MatrixXd diff_occ(p, n_occ);
  for (int j = 0; j < n_occ; ++j) {
    diff_occ.col(j) = occ_eig.row(j).transpose() - mu_eig;
  }
  
  Eigen::MatrixXd diff_m(p, n_m);
  for (int j = 0; j < n_m; ++j) {
    diff_m.col(j) = m_eig.row(j).transpose() - mu_eig;
  }
  
  // Solve triangular systems for all columns at once
  Eigen::MatrixXd y_occ = L_eig.triangularView<Eigen::Lower>().solve(diff_occ);
  Eigen::MatrixXd y_m   = L_eig.triangularView<Eigen::Lower>().solve(diff_m);
  
  // Sum of squared column norms
  double sum_q1 = y_occ.colwise().squaredNorm().sum();
  
  // a_j = -0.5 * ||y_j||^2  (as array, 1 x n_m)
  Eigen::ArrayXd a = -0.5 * y_m.colwise().squaredNorm().array();
  // Stable log-sum-exp
  double max_a = a.maxCoeff();
  double sum_exp = (a - max_a).exp().sum();
  double log_sum_exp = max_a + std::log(sum_exp);
  
  double neg_log = 0.5 * sum_q1 + static_cast<double>(n_occ) * log_sum_exp;
  
  if (!std::isfinite(neg_log)) {
    neg_log = nicher::OPTIM_PENALTY;
  }
  
  return neg_log;
}