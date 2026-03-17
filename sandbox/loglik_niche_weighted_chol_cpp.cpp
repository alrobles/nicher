#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <limits>

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

// [[Rcpp::export]]
double loglik_niche_weighted_chol_cpp(NumericVector mu,
                                      NumericMatrix L,          // factor Cholesky triangular inferior de Sigma
                                      NumericMatrix env_occ,    // puntos de presencia (n_occ x p)
                                      NumericMatrix env_den,    // puntos de M para denominador (n_den x p)
                                      NumericVector w_occ,      // pesos KDE en presencias (n_occ)
                                      NumericVector w_den) {    // pesos KDE en puntos de denominador (n_den)
  int p = L.nrow();
  int n_occ = env_occ.nrow();
  int n_den = env_den.nrow();
  
  // Validaciones
  if (L.ncol() != p) stop("L must be square");
  if (env_occ.ncol() != p) stop("env_occ must have p columns");
  if (env_den.ncol() != p) stop("env_den must have p columns");
  if (mu.size() != p) stop("mu must have length p");
  if (n_occ <= 0 || n_den <= 0) stop("n_occ and n_den must be >0");
  if (w_occ.size() != n_occ) stop("w_occ length must match n_occ");
  if (w_den.size() != n_den) stop("w_den length must match n_den");
  
  // Mapear a Eigen (sin copia)
  Eigen::Map<Eigen::VectorXd> mu_eig(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(mu));
  Eigen::Map<Eigen::MatrixXd> L_eig(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(L));
  Eigen::Map<Eigen::MatrixXd> occ_eig(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(env_occ));
  Eigen::Map<Eigen::MatrixXd> den_eig(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(env_den));
  
  // Map pesos (los trataremos como vectores de Eigen para operaciones array)
  Eigen::Map<Eigen::VectorXd> w_occ_eig(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(w_occ));
  Eigen::Map<Eigen::VectorXd> w_den_eig(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(w_den));
  
  // --- Calcular diferencias y resolver sistemas triangulares ---
  // Matriz de diferencias para presencias: p x n_occ
  Eigen::MatrixXd diff_occ(p, n_occ);
  for (int j = 0; j < n_occ; ++j) {
    diff_occ.col(j) = occ_eig.row(j).transpose() - mu_eig;
  }
  
  // Matriz de diferencias para denominador: p x n_den
  Eigen::MatrixXd diff_den(p, n_den);
  for (int j = 0; j < n_den; ++j) {
    diff_den.col(j) = den_eig.row(j).transpose() - mu_eig;
  }
  
  // Resolver L * Y = diff para todas las columnas
  Eigen::MatrixXd y_occ = L_eig.triangularView<Eigen::Lower>().solve(diff_occ);
  Eigen::MatrixXd y_den = L_eig.triangularView<Eigen::Lower>().solve(diff_den);
  
  // Suma de Mahalanobis para presencias
  double sum_q1 = y_occ.colwise().squaredNorm().sum();
  
  // q2 para denominador (como array)
  Eigen::ArrayXd q2 = y_den.colwise().squaredNorm().array();
  
  // Log-pesos (tomamos logaritmo natural)
  Eigen::ArrayXd log_w_occ = w_occ_eig.array().log();
  Eigen::ArrayXd log_w_den = w_den_eig.array().log();
  
  // --- Término log-sum-exp para el denominador ---
  // a_j = log(w_den_j) - 0.5 * q2_j
  Eigen::ArrayXd a = log_w_den - 0.5 * q2;
  
  double max_a = a.maxCoeff();
  double sum_exp = (a - max_a).exp().sum();
  double log_sum_exp = max_a + std::log(sum_exp);
  
  // --- Log-verosimilitud negativa ---
  double neg_log = 0.5 * sum_q1 - log_w_occ.sum() + static_cast<double>(n_occ) * log_sum_exp;
  
  return neg_log;
}