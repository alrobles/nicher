#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector kde_gaussian_2d_cpp(const NumericMatrix& x, const NumericMatrix& data) {
  int n_eval = x.nrow();
  int n_data = data.nrow();
  int p = data.ncol();
  
  if (p != 2) stop("kde_gaussian_2d_cpp only works for p=2");
  if (n_data <= 0) stop("data must have n>0");
  if (x.ncol() != p) stop("x and data must have the same number of columns");
  
  // Calcular medias y varianzas
  double mean1 = 0.0, mean2 = 0.0;
  for (int i = 0; i < n_data; ++i) {
    mean1 += data(i, 0);
    mean2 += data(i, 1);
  }
  mean1 /= n_data;
  mean2 /= n_data;
  
  double var1 = 0.0, var2 = 0.0;
  for (int i = 0; i < n_data; ++i) {
    double d1 = data(i, 0) - mean1;
    double d2 = data(i, 1) - mean2;
    var1 += d1 * d1;
    var2 += d2 * d2;
  }
  var1 /= n_data;
  var2 /= n_data;
  
  double s1 = std::sqrt(var1);
  double s2 = std::sqrt(var2);
  if (s1 <= 0) s1 = 1e-8;
  if (s2 <= 0) s2 = 1e-8;
  
  double nf = std::pow(static_cast<double>(n_data), -1.0 / 6.0); // p+4 = 6
  double bw1 = s1 * nf;
  double bw2 = s2 * nf;
  double inv_bw1_2 = 1.0 / (bw1 * bw1);
  double inv_bw2_2 = 1.0 / (bw2 * bw2);
  double prod_bw = bw1 * bw2;
  if (prod_bw <= 0) prod_bw = 1e-32;
  
  double log2pi = std::log(2.0 * M_PI);
  double log_const = -0.5 * 2 * log2pi - std::log(prod_bw);
  double const_kernel = std::exp(log_const);
  
  NumericVector dens(n_eval);
  for (int i = 0; i < n_eval; ++i) {
    double xi1 = x(i, 0);
    double xi2 = x(i, 1);
    double sum = 0.0;
    for (int j = 0; j < n_data; ++j) {
      double d1 = data(j, 0) - xi1;
      double d2 = data(j, 1) - xi2;
      double q = d1 * d1 * inv_bw1_2 + d2 * d2 * inv_bw2_2;
      sum += std::exp(-0.5 * q);
    }
    dens[i] = const_kernel * (sum / n_data);
  }
  return dens;
}