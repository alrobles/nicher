// src/loglik_niche_presence_only.cpp
#include "nicher_types.h"

using namespace Rcpp;

// [[Rcpp::export]]
double loglik_niche_presence_only_cpp(NumericVector mu,
                                      NumericMatrix L,
                                      NumericMatrix env_occ) {
  int p = L.nrow();
  int n_occ = env_occ.nrow();
  
  // Basic input validation
  if (L.ncol() != p) stop("L must be square");
  if (env_occ.ncol() != p) stop("env_occ must have p columns");
  if (mu.size() != p) stop("mu must have length p");
  if (n_occ <= 0) stop("n_occ must be >0");
  
  // Map to Eigen
  Eigen::Map<Eigen::VectorXd> mu_eig(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(mu));
  Eigen::Map<Eigen::MatrixXd> L_eig(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(L));
  Eigen::Map<Eigen::MatrixXd> occ_eig(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(env_occ));
  
  // Compute sum of Mahalanobis distances
  double sum_q = nicher::sum_mahalanobis_sq(occ_eig, mu_eig, L_eig);
  
  // Log-determinant: 2 * sum(log(diag(L))) because |Sigma| = (det(L))^2 and det(L) = prod(diag(L))
  double log_det = 0.0;
  for (int i = 0; i < p; ++i) {
    log_det += std::log(L_eig(i, i));
  }
  log_det *= 2.0;  // because |Sigma| = (prod diag(L))^2, so log|Sigma| = 2 * sum(log diag(L))
  
  // Negative log-likelihood (excluding constants that do not depend on parameters,
  // but including the normalisation term):
  // Full likelihood: (2pi)^(-n*p/2) * |Sigma|^(-n/2) * exp(-0.5 * sum_q)
  // Log and negate: 0.5 * n * (p*log(2pi) + log|Sigma|) + 0.5 * sum_q
  // For minimization we can drop constants (p*log(2pi)) but keep log|Sigma|.
  double n = static_cast<double>(n_occ);
  double neg_log = 0.5 * n * log_det + 0.5 * sum_q;
  
  if (!std::isfinite(neg_log)) {
    neg_log = nicher::OPTIM_PENALTY;
  }
  
  return neg_log;
}