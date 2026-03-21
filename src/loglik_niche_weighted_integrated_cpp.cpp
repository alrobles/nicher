// src/loglik_niche_weighted_integrated.cpp
#include "nicher_types.h"

using namespace Rcpp;
using namespace nicher; // para usar nuestras funciones KDE

// [[Rcpp::export]]
double loglik_niche_weighted_integrated_cpp(
    NumericVector mu,
    NumericMatrix L,
    NumericMatrix env_occ,
    NumericMatrix env_m,
    Nullable<IntegerVector> den_idx = R_NilValue,
    Nullable<IntegerVector> kde_idx = R_NilValue,
    Nullable<NumericVector> precomp_w_den = R_NilValue,   // nuevo: pesos precalculados para el denominador
    bool neg = true) {
  
  int p = L.nrow();
  int n_occ = env_occ.nrow();
  int n_m = env_m.nrow();
  
  // Validaciones
  if (L.ncol() != p) stop("L must be square");
  if (env_occ.ncol() != p) stop("env_occ must have p columns");
  if (env_m.ncol() != p) stop("env_m must have p columns");
  if (mu.size() != p) stop("mu must have length p");
  if (n_occ <= 0 || n_m <= 0) stop("n_occ and n_m must be >0");
  
  // Mapear a Eigen
  Eigen::Map<Eigen::VectorXd> mu_eig(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(mu));
  Eigen::Map<Eigen::MatrixXd> L_eig(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(L));
  Eigen::Map<Eigen::MatrixXd> occ_eig(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(env_occ));
  Eigen::Map<Eigen::MatrixXd> m_eig(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(env_m));
  
  // --- Construir subconjuntos según índices ---
  Eigen::MatrixXd M_den, M_kde;
  
  // Denominador
  if (den_idx.isNotNull()) {
    IntegerVector idx(den_idx);
    int n_den = idx.size();
    M_den.resize(n_den, p);
    for (int i = 0; i < n_den; ++i) {
      int row = idx[i] - 1; // R tiene índices 1-based
      if (row < 0 || row >= n_m) stop("den_idx out of range");
      M_den.row(i) = m_eig.row(row);
    }
  } else {
    M_den = m_eig; // Nota: esto copia, pero es inevitable para tener un objeto Eigen propio
  }
  
  // KDE reference
  if (kde_idx.isNotNull()) {
    IntegerVector idx(kde_idx);
    int n_kde = idx.size();
    M_kde.resize(n_kde, p);
    for (int i = 0; i < n_kde; ++i) {
      int row = idx[i] - 1;
      if (row < 0 || row >= n_m) stop("kde_idx out of range");
      M_kde.row(i) = m_eig.row(row);
    }
  } else {
    M_kde = m_eig;
  }
  
  // --- Calcular pesos KDE ---
  Eigen::VectorXd w_occ, w_den;
  
  // Pesos para las presencias (siempre se calculan)
  if (p == 2) {
    w_occ = nicher::kde_2d(occ_eig, M_kde);
  } else {
    w_occ = nicher::kde_eigen(occ_eig, M_kde);
  }
  
  // Pesos para el denominador: usar precomputado si se proporciona, sino calcular
  if (precomp_w_den.isNotNull()) {
    // Convertir el vector Numeric a Eigen::VectorXd
    w_den = Rcpp::as<Eigen::Map<Eigen::VectorXd> >(precomp_w_den);
    // Verificar que la longitud coincide con el número de filas de M_den
    if (w_den.size() != M_den.rows()) {
      stop("Length of precomp_w_den (%d) does not match number of denominator rows (%d)",
           w_den.size(), M_den.rows());
    }
  } else {
    // Calcular pesos normalmente
    if (p == 2) {
      w_den = nicher::kde_2d(M_den, M_kde);
    } else {
      w_den = nicher::kde_eigen(M_den, M_kde);
    }
  }
  
  // --- Distancias de Mahalanobis usando factor L ---
  // Matriz de diferencias para presencias
  Eigen::MatrixXd diff_occ(p, n_occ);
  for (int j = 0; j < n_occ; ++j) {
    diff_occ.col(j) = occ_eig.row(j).transpose() - mu_eig;
  }
  Eigen::MatrixXd y_occ = L_eig.triangularView<Eigen::Lower>().solve(diff_occ);
  double sum_q1 = y_occ.colwise().squaredNorm().sum();
  
  // Matriz de diferencias para denominador
  int n_den = M_den.rows();
  Eigen::MatrixXd diff_den(p, n_den);
  for (int j = 0; j < n_den; ++j) {
    diff_den.col(j) = M_den.row(j).transpose() - mu_eig;
  }
  Eigen::MatrixXd y_den = L_eig.triangularView<Eigen::Lower>().solve(diff_den);
  Eigen::ArrayXd q2 = y_den.colwise().squaredNorm().array();
  
  // --- Log-verosimilitud ---
  // Clamp KDE weights to a minimum positive value to prevent log(0) = -Inf
  // and subsequent NaN propagation through the optimizer
  const double min_kde_weight = 1e-300;
  Eigen::ArrayXd log_w_occ = w_occ.array().max(min_kde_weight).log();
  Eigen::ArrayXd log_w_den = w_den.array().max(min_kde_weight).log();
  Eigen::ArrayXd a = log_w_den - 0.5 * q2;
  
  double max_a = a.maxCoeff();
  double sum_exp = (a - max_a).exp().sum();
  double log_sum_exp = max_a + std::log(sum_exp);
  
  double neg_log = 0.5 * sum_q1 - log_w_occ.sum() + static_cast<double>(n_occ) * log_sum_exp;
  
  // Guard against NaN/Inf: return a large finite penalty so the optimizer
  // can recover instead of crashing with a type conversion error
  if (!std::isfinite(neg_log)) {
    neg_log = 1e300;
  }
  
  return neg ? neg_log : -neg_log;
}