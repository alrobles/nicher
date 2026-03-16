// src/loglik_niche_cpp.cpp
#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
double loglik_niche_cpp(NumericVector mu,
                        NumericMatrix S,          // covariance matrix (p x p)
                        NumericMatrix env_occ,    // n_occ x p
                        NumericMatrix env_m) {    // n_m   x p
  const int p = S.nrow();
  const int n_occ = env_occ.nrow();
  const int n_m = env_m.nrow();
  
  // Validaciones básicas
  if (S.ncol() != p) stop("S must be square p x p");
  if (env_occ.ncol() != p) stop("env_occ must have p columns");
  if (env_m.ncol() != p) stop("env_m must have p columns");
  if (mu.size() != p) stop("mu must have length p");
  if (n_occ <= 0 || n_m <= 0) stop("env_occ/env_m must have n>0");
  
  // Mapear objetos de R a Eigen (sin copia)
  Map<MatrixXd> S_eig(as<Map<MatrixXd> >(S));           // p x p
  Map<VectorXd> mu_eig(as<Map<VectorXd> >(mu));         // p
  Map<MatrixXd> occ_eig(as<Map<MatrixXd> >(env_occ));   // n_occ x p
  Map<MatrixXd> m_eig(as<Map<MatrixXd> >(env_m));       // n_m x p
  
  // Descomposición de Cholesky de S
  LLT<MatrixXd> llt(S_eig);
  if (llt.info() != Success) {
    stop("Cholesky decomposition failed: covariance matrix not positive definite");
  }
  
  // Extraer el factor triangular inferior L (S = L L^T)
  MatrixXd L = llt.matrixL();
  
  // Función lambda para calcular la distancia de Mahalanobis al cuadrado
  auto mahalanobis_sq = [&](const VectorXd& x) -> double {
    VectorXd d = x - mu_eig;               // diferencia
    // Resolver L * y = d  (sustitución hacia adelante)
    VectorXd y = L.triangularView<Lower>().solve(d);
    return y.squaredNorm();                 // y^T y
  };
  
  // Suma de D^2 para ocurrencias
  double sum_q1 = 0.0;
  for (int i = 0; i < n_occ; ++i) {
    sum_q1 += mahalanobis_sq(occ_eig.row(i).transpose());
  }
  
  // Calcular a_j = -0.5 * D^2 para M y hacer log-sum-exp estable
  std::vector<double> a(n_m);
  double max_a = -std::numeric_limits<double>::infinity();
  for (int j = 0; j < n_m; ++j) {
    double d2 = mahalanobis_sq(m_eig.row(j).transpose());
    double aj = -0.5 * d2;
    a[j] = aj;
    if (aj > max_a) max_a = aj;
  }
  
  double sum_exp = 0.0;
  for (int j = 0; j < n_m; ++j) {
    sum_exp += std::exp(a[j] - max_a);
  }
  double log_sum_exp = max_a + std::log(sum_exp);
  
  // -log L (sin constantes que se cancelan)
  double neg_log = 0.5 * sum_q1 + static_cast<double>(n_occ) * log_sum_exp;
  return neg_log;
}