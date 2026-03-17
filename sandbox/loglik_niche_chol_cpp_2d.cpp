#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double loglik_niche_chol_cpp_2d(NumericVector mu,
                                NumericMatrix L,          // factor Cholesky 2x2 inferior
                                NumericMatrix env_occ,
                                NumericMatrix env_m) {
  // Validaciones
  if (L.nrow() != 2 || L.ncol() != 2) stop("L must be 2x2");
  if (env_occ.ncol() != 2) stop("env_occ must have 2 columns");
  if (env_m.ncol() != 2) stop("env_m must have 2 columns");
  if (mu.size() != 2) stop("mu must have length 2");
  int n_occ = env_occ.nrow();
  int n_m = env_m.nrow();
  if (n_occ <= 0 || n_m <= 0) stop("env_occ/env_m must have n>0");
  
  // Extraer elementos de L (triangular inferior)
  double L11 = L(0,0);
  double L21 = L(1,0);
  double L22 = L(1,1);
  double mu1 = mu[0];
  double mu2 = mu[1];
  
  // Precomputar inversos de diagonales para acelerar (opcional)
  double invL11 = 1.0 / L11;
  double invL22 = 1.0 / L22;
  
  // Suma de Mahalanobis para ocurrencias
  double sum_q1 = 0.0;
  for (int i = 0; i < n_occ; ++i) {
    double d1 = env_occ(i,0) - mu1;
    double d2 = env_occ(i,1) - mu2;
    // Resolver L y = d  =>  y1 = d1 / L11,   y2 = (d2 - L21*y1) / L22
    double y1 = d1 * invL11;
    double y2 = (d2 - L21 * y1) * invL22;
    sum_q1 += y1*y1 + y2*y2;
  }
  
  // Vector para almacenar a_j = -0.5 * D^2
  std::vector<double> a(n_m);
  double max_a = -std::numeric_limits<double>::infinity();
  
  for (int j = 0; j < n_m; ++j) {
    double d1 = env_m(j,0) - mu1;
    double d2 = env_m(j,1) - mu2;
    double y1 = d1 * invL11;
    double y2 = (d2 - L21 * y1) * invL22;
    double d2_val = y1*y1 + y2*y2;
    double aj = -0.5 * d2_val;
    a[j] = aj;
    if (aj > max_a) max_a = aj;
  }
  
  // Log-sum-exp estable
  double sum_exp = 0.0;
  for (int j = 0; j < n_m; ++j) {
    sum_exp += std::exp(a[j] - max_a);
  }
  double log_sum_exp = max_a + std::log(sum_exp);
  
  double neg_log = 0.5 * sum_q1 + static_cast<double>(n_occ) * log_sum_exp;
  return neg_log;
}