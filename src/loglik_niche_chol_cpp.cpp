#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

// [[Rcpp::export]]
double loglik_niche_chol_cpp(NumericVector mu,
                             NumericMatrix L,          // factor Cholesky triangular inferior de Sigma
                             NumericMatrix env_occ,
                             NumericMatrix env_m) {
  int p = L.nrow();
  int n_occ = env_occ.nrow();
  int n_m = env_m.nrow();
  
  // Validaciones
  if (L.ncol() != p) stop("L must be square");
  if (env_occ.ncol() != p) stop("env_occ must have p columns");
  if (env_m.ncol() != p) stop("env_m must have p columns");
  if (mu.size() != p) stop("mu must have length p");
  if (n_occ <= 0 || n_m <= 0) stop("env_occ/env_m must have n>0");
  
  // Mapear a Eigen (sin copia)
  Eigen::Map<Eigen::VectorXd> mu_eig(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(mu));
  Eigen::Map<Eigen::MatrixXd> L_eig(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(L));
  Eigen::Map<Eigen::MatrixXd> occ_eig(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(env_occ));
  Eigen::Map<Eigen::MatrixXd> m_eig(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(env_m));
  
  // Función para distancia de Mahalanobis al cuadrado usando L
  auto mahalanobis_sq = [&](const Eigen::VectorXd& x) -> double {
    Eigen::VectorXd d = x - mu_eig;
    // Resolver L * y = d (sustitución hacia adelante)
    Eigen::VectorXd y = L_eig.triangularView<Eigen::Lower>().solve(d);
    return y.squaredNorm();
  };
  
  double sum_q1 = 0.0;
  for (int i = 0; i < n_occ; ++i) {
    sum_q1 += mahalanobis_sq(occ_eig.row(i).transpose());
  }
  
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
  
  double neg_log = 0.5 * sum_q1 + static_cast<double>(n_occ) * log_sum_exp;
  return neg_log;
}