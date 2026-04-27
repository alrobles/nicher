// niche_suitability_cpp.cpp — RcppParallel kernel for habitat-suitability
// evaluation. Single-pass over a column-major (n_loc x p) flat buffer; per
// pixel computes the standardized Gaussian
//
//   S(x) = exp( -1/2 * (x - mu)^T Sigma^{-1} (x - mu) )
//
// (Jimenez et al. 2022, Eq. 2 + standardization paragraph). The Cholesky
// decomposition Sigma = L L^T and its inverse L_inv are precomputed R-side
// once, so the kernel only does a triangular-solve per pixel (~p^2/2 flops).
//
// Buffer layout matches xsdm-devel recipe 03 (terra block-wise I/O):
//
//   env_vec[l + n_loc * k] = value of variable k at location l
//
// References:
//   - xsdm-devel docs/methodology/03-terra-blockwise-io.md
//   - xsdm-devel docs/methodology/04-ucminfcpp-xptr.md (try/catch wrapper)
//   - Jimenez et al. (2022) Ecological Modelling 464, 109823.

// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>

#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;

namespace {

// Per-pixel worker. Reads from a flat const double* buffer with no R API
// calls in the parallel section (RcppParallel safety contract).
struct NicheSuitabilityWorker : public Worker {
  const double* env_ptr;     // length n_loc * p, column-major
  int           n_loc;
  int           p;
  const double* mu_ptr;      // length p
  const double* Linv_ptr;    // p * p, column-major lower-triangular
  bool          return_log;
  double*       out_ptr;     // length n_loc

  NicheSuitabilityWorker(const double* env_ptr_, int n_loc_, int p_,
                         const double* mu_ptr_, const double* Linv_ptr_,
                         bool return_log_, double* out_ptr_)
    : env_ptr(env_ptr_), n_loc(n_loc_), p(p_),
      mu_ptr(mu_ptr_), Linv_ptr(Linv_ptr_),
      return_log(return_log_), out_ptr(out_ptr_) {}

  void operator()(std::size_t begin, std::size_t end) {
    // Stack-allocated scratch — small (p typically <= 6).
    constexpr int MAX_P_STACK = 32;
    double r_stack[MAX_P_STACK];
    double z_stack[MAX_P_STACK];
    std::vector<double> r_heap, z_heap;
    double* r;
    double* z;
    if (p <= MAX_P_STACK) {
      r = r_stack;
      z = z_stack;
    } else {
      r_heap.assign(p, 0.0);
      z_heap.assign(p, 0.0);
      r = r_heap.data();
      z = z_heap.data();
    }

    for (std::size_t l = begin; l < end; ++l) {
      // Wrap each pixel in try/catch (xsdm-devel recipe 04 step 8) so
      // pathological inputs surface as NA rather than aborting the
      // worker thread.
      try {
        // r = x - mu
        bool any_nan = false;
        for (int i = 0; i < p; ++i) {
          double xi = env_ptr[static_cast<std::size_t>(l) +
                              static_cast<std::size_t>(n_loc) * i];
          if (!R_FINITE(xi)) { any_nan = true; break; }
          r[i] = xi - mu_ptr[i];
        }
        if (any_nan) {
          out_ptr[l] = NA_REAL;
          continue;
        }

        // z = L^{-1} r  (lower-triangular forward solve, expanded form)
        // L_inv stored column-major: Linv_ptr[i + p*j] = L_inv(i, j),
        // upper triangle is zero.
        double m = 0.0;
        for (int i = 0; i < p; ++i) {
          double zi = 0.0;
          for (int j = 0; j <= i; ++j) {
            zi += Linv_ptr[i + static_cast<std::size_t>(p) * j] * r[j];
          }
          z[i] = zi;
          m   += zi * zi;
        }

        const double v = -0.5 * m;
        out_ptr[l] = return_log ? v : std::exp(v);
      } catch (...) {
        out_ptr[l] = NA_REAL;
      }
    }
  }
};

}  // namespace

//' Habitat-suitability kernel (parallel, zero-copy)
//'
//' Evaluates the standardized multivariate-normal density
//' \eqn{S(x) = \exp(-\tfrac12 (x-\mu)^\top \Sigma^{-1} (x-\mu))} for each
//' row of a column-major flat environmental buffer. Uses
//' \pkg{RcppParallel}'s \code{parallelFor} over locations.
//'
//' @param env_dat_vec Numeric vector. Flat column-major buffer of length
//'   \code{n_loc * p}; entry for location \code{l} variable \code{k} sits
//'   at index \code{l + n_loc * k}.
//' @param env_dat_dims Integer vector \code{c(n_loc, p)}.
//' @param mu Numeric vector of length \code{p}: niche centroid.
//' @param L_inv Numeric matrix \code{p x p}: inverse of the lower
//'   Cholesky factor of \eqn{\Sigma}. Precomputed R-side once.
//' @param return_log Logical. If \code{FALSE} (default) returns
//'   suitability in \eqn{(0, 1]}; if \code{TRUE} returns
//'   \eqn{\log S(x) \le 0}.
//' @param num_threads Integer. \code{0} (default) leaves
//'   \code{RcppParallel}'s global thread state unchanged.
//'
//' @return Numeric vector of length \code{n_loc}.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector niche_suitability_cpp(
    Rcpp::NumericVector env_dat_vec,
    Rcpp::IntegerVector env_dat_dims,
    Rcpp::NumericVector mu,
    Rcpp::NumericMatrix L_inv,
    bool return_log = false,
    int  num_threads = 0
) {
  if (env_dat_dims.size() != 2) {
    Rcpp::stop("env_dat_dims must be c(n_loc, p)");
  }
  const int n_loc = env_dat_dims[0];
  const int p     = env_dat_dims[1];
  if (n_loc < 0 || p < 0) {
    Rcpp::stop("env_dat_dims entries must be non-negative");
  }
  if (env_dat_vec.size() != static_cast<R_xlen_t>(n_loc) * p) {
    Rcpp::stop("env_dat_vec length must equal n_loc * p");
  }
  if (mu.size() != p) {
    Rcpp::stop("mu must have length p");
  }
  if (L_inv.nrow() != p || L_inv.ncol() != p) {
    Rcpp::stop("L_inv must be a p x p matrix");
  }

  Rcpp::NumericVector out(n_loc);
  if (n_loc == 0) return out;

  // Thread management — same pattern as xsdm-devel src/log_prob_detect.cpp:
  // route through RcppParallel's R namespace so we don't depend on the
  // private TBB headers. num_threads = 0 → leave the global state alone.
  int  old_threads      = -1;
  bool restore_threads  = false;
  Rcpp::Environment rcppPar  = Rcpp::Environment::namespace_env("RcppParallel");
  Rcpp::Function    defaultNT = rcppPar["defaultNumThreads"];
  Rcpp::Function    setTO     = rcppPar["setThreadOptions"];
  if (num_threads > 0) {
    old_threads     = Rcpp::as<int>(defaultNT());
    restore_threads = true;
    setTO(Rcpp::Named("numThreads") = num_threads);
  }

  NicheSuitabilityWorker w(
      REAL(env_dat_vec), n_loc, p,
      REAL(mu), REAL(L_inv),
      return_log, REAL(out));

  // Grain size keeps small rasters from being chopped into single-pixel
  // tasks (parallel overhead would dominate).
  RcppParallel::parallelFor(0, static_cast<std::size_t>(n_loc), w, 1024);

  if (restore_threads) {
    setTO(Rcpp::Named("numThreads") = old_threads);
  }

  return out;
}
