#include <Rcpp.h>
using namespace Rcpp;

// Calcula sd por columna (sin NAs)
static void col_sd(const NumericMatrix& X, NumericVector& sd) {
  const int n = X.nrow();
  const int p = X.ncol();
  for (int j = 0; j < p; ++j) {
    double s1 = 0.0, s2 = 0.0;
    for (int i = 0; i < n; ++i) {
      const double v = X(i, j);
      s1 += v;
      s2 += v * v;
    }
    const double m = s1 / n;
    double var = (s2 / n) - m * m;
    if (var < 0.0) var = 0.0;
    sd[j] = std::sqrt(var);
    if (!R_finite(sd[j]) || sd[j] <= 0.0) sd[j] = 1e-8; // evita sd=0
  }
}

// [[Rcpp::export]]
NumericVector kde_gaussian_rcpp(const NumericMatrix& x,
                                      const NumericMatrix& data) {
  const int n_eval = x.nrow();
  const int n_data = data.nrow();
  const int p = data.ncol();
  if (n_data <= 0 || p <= 0) stop("data must have n>0 and p>0");
  if (x.ncol() != p) stop("x and data must have the same number of columns (p)");
  
  // Scott diagonal: h_j = sd_j * n^(-1/(p+4))
  NumericVector s(p);
  col_sd(data, s);
  const double expo = -1.0 / (p + 4.0);
  const double nf   = std::pow((double)n_data, expo);
  
  NumericVector bw(p), inv_bwsq(p);
  double prod_bw = 1.0;
  for (int j = 0; j < p; ++j) {
    bw[j] = s[j] * nf;
    if (!R_finite(bw[j]) || bw[j] <= 0.0) bw[j] = 1e-8;
    inv_bwsq[j] = 1.0 / (bw[j] * bw[j]);
    prod_bw *= bw[j];
  }
  if (!R_finite(prod_bw) || prod_bw <= 0.0) prod_bw = 1e-32;
  
  // Constante del kernel gaussiano con H = diag(bw^2)
  const double log2pi = std::log(2.0 * M_PI);
  const double log_const = -0.5 * p * log2pi - std::log(prod_bw);
  const double const_kernel = std::exp(log_const);
  
  NumericVector dens(n_eval);
  for (int i = 0; i < n_eval; ++i) {
    double acc = 0.0;
    for (int r = 0; r < n_data; ++r) {
      double q = 0.0;
      for (int j = 0; j < p; ++j) {
        const double dx = data(r, j) - x(i, j);
        q += dx * dx * inv_bwsq[j];
      }
      acc += std::exp(-0.5 * q);
    }
    dens[i] = const_kernel * (acc / (double)n_data);
  }
  return dens;
}