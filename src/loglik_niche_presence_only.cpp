// src/loglik_niche_presence_only.cpp
#include "nicher_types.h"

using namespace Rcpp;  // <-- AÑADIR ESTA LÍNEA

// [[Rcpp::export]]
double loglik_niche_presence_only_cpp(NumericVector mu,
                                      NumericMatrix L,
                                      NumericMatrix env_occ) {
  int p = L.nrow();
  int n_occ = env_occ.nrow();
  
  // Validaciones básicas
  if (L.ncol() != p) stop("L must be square");
  if (env_occ.ncol() != p) stop("env_occ must have p columns");
  if (mu.size() != p) stop("mu must have length p");
  if (n_occ <= 0) stop("n_occ must be >0");
  
  // Mapear a Eigen
  Eigen::Map<Eigen::VectorXd> mu_eig(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(mu));
  Eigen::Map<Eigen::MatrixXd> L_eig(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(L));
  Eigen::Map<Eigen::MatrixXd> occ_eig(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(env_occ));
  
  // Calcular suma de Mahalanobis
  double sum_q = nicher::sum_mahalanobis_sq(occ_eig, mu_eig, L_eig);
  
  // Log-determinante: 2 * suma(log(diag(L))) porque |Sigma| = (det(L))^2 y det(L) = prod(diag(L))
  double log_det = 0.0;
  for (int i = 0; i < p; ++i) {
    log_det += std::log(L_eig(i, i));
  }
  log_det *= 2.0;  // porque |Sigma| = (prod diag(L))^2, entonces log|Sigma| = 2 * sum(log diag(L))
  
  // Log-verosimilitud negativa (sin constantes que no dependen de parámetros, pero incluyendo el término de normalización)
  // La verosimilitud completa: (2pi)^(-n*p/2) * |Sigma|^(-n/2) * exp(-0.5 * sum_q)
  // Tomando log y negativo: 0.5 * n * (p*log(2pi) + log|Sigma|) + 0.5 * sum_q
  // Para minimización, podemos omitir constantes (p*log(2pi)) pero mantener log|Sigma|.
  double n = static_cast<double>(n_occ);
  double neg_log = 0.5 * n * log_det + 0.5 * sum_q;
  
  if (!std::isfinite(neg_log)) {
    neg_log = nicher::OPTIM_PENALTY;
  }
  
  return neg_log;
}