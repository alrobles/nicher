// src/niche_obj.cpp
//
// Builds the `Rcpp::XPtr<ucminf::ObjFun>` consumed by `ucminfcpp::ucminf_xptr`.
//
// Compared with the previous version, the closure now:
//   * routes through the pure-C++ math-scale kernels in
//     `loglik_niche_math_cpp.cpp` (no Rcpp::NumericMatrix copies inside the
//     optimizer's inner loop),
//   * supports `grad = "analytic"` for the weighted likelihood (hybrid
//     analytic over mu / log_sigma + FD over the C-vine v block), in
//     addition to the previous central / forward FD modes,
//   * accepts precomputed KDE weights for BOTH the presence side
//     (`precomp_w_occ`) and the M-denominator side (`precomp_w_den`), so
//     that the KDE computation is hoisted out of the optimization loop, and
//   * wraps the entire closure in try/catch boundaries so a bad evaluation
//     surfaces as `f = +Inf, g = 0` rather than aborting R inside
//     `ucminf::minimize_direct<F>`. See `xsdm-devel/docs/methodology/
//     04-ucminfcpp-xptr.md` step 8 for the rationale.

#include "nicher_types.h"
#include <functional>
#include <vector>
#include <cmath>
#include <memory>
#include <string>
#include <stdexcept>
using namespace Rcpp;

// Forward declaration of the existing (mu, L)-scale unweighted kernel; we
// keep the unweighted path for backward compatibility with existing tests.
double loglik_niche_chol_cpp(
    NumericVector mu, NumericMatrix L,
    NumericMatrix env_occ, NumericMatrix env_m);

NumericMatrix cvine_cholesky(NumericVector v, int d, double eta);

// ---------------------------------------------------------------------------
// Closure state
// ---------------------------------------------------------------------------

namespace {

enum LikType { LIK_UNWEIGHTED = 0, LIK_WEIGHTED = 1, LIK_PRESENCE_ONLY = 2 };
enum GradMode { GRAD_ANALYTIC = 0, GRAD_CENTRAL = 1, GRAD_FORWARD = 2 };

struct XptrClosureState {
  int p;
  double eta;
  LikType  lik_type;
  GradMode grad_mode;
  double gradstep_rel;
  double gradstep_abs;

  // Snapshot owned at construction. Per-call work uses these directly --
  // no allocations on the optimizer hot path.
  Eigen::MatrixXd env_occ;          // n_occ x p
  Eigen::MatrixXd env_m_full;       // for unweighted only (n_m x p)
  Eigen::MatrixXd M_den;             // weighted: denominator slice (n_den x p)
  Eigen::VectorXd w_occ;             // weighted: KDE at presence rows
  Eigen::VectorXd w_den;             // weighted: KDE at denominator rows
};

using NicheObjFunType = std::function<void(const std::vector<double>&,
                                           std::vector<double>&,
                                           double&)>;

// Helper for the unweighted (legacy) path: build the (mu, L) representation
// and call the existing kernel. Kept verbatim from the prior implementation
// so that benchmark tests against `loglik_niche_chol_cpp` stay green.
static double eval_unweighted_legacy(const std::vector<double>& x,
                                     const XptrClosureState& s) {
  const int p = s.p;
  const int n_v = p * (p - 1) / 2;
  if ((int)x.size() != 2 * p + n_v) Rcpp::stop("theta length mismatch");

  NumericVector mu(p), v(n_v);
  std::vector<double> sigma(p);
  for (int i = 0; i < p; ++i) {
    mu[i] = x[i];
    sigma[i] = std::exp(x[p + i]);
  }
  for (int k = 0; k < n_v; ++k) v[k] = x[2 * p + k];

  NumericMatrix L_corr = cvine_cholesky(v, p, s.eta);
  NumericMatrix L_cov(p, p);
  for (int j = 0; j < p; ++j)
    for (int k = 0; k <= j; ++k)
      L_cov(j, k) = sigma[j] * L_corr(j, k);

  NumericMatrix occ_m(s.env_occ.rows(), s.env_occ.cols());
  std::memcpy(occ_m.begin(), s.env_occ.data(),
              sizeof(double) * s.env_occ.size());
  NumericMatrix m_m(s.env_m_full.rows(), s.env_m_full.cols());
  std::memcpy(m_m.begin(), s.env_m_full.data(),
              sizeof(double) * s.env_m_full.size());

  return loglik_niche_chol_cpp(mu, L_cov, occ_m, m_m);
}

// Compute the objective value (no gradient) for any likelihood type.
static double eval_value(const std::vector<double>& x,
                         const XptrClosureState& s) {
  switch (s.lik_type) {
    case LIK_PRESENCE_ONLY:
      return nicher::loglik_niche_math_presence_only_eigen(
          x.data(), (int)x.size(), s.env_occ, s.eta);
    case LIK_WEIGHTED:
      return nicher::loglik_niche_math_weighted_eigen(
          x.data(), (int)x.size(), s.env_occ, s.M_den,
          s.w_occ, s.w_den, s.eta);
    case LIK_UNWEIGHTED:
    default:
      return eval_unweighted_legacy(x, s);
  }
}

// Central / forward finite-difference gradient (used for presence_only and
// unweighted, and as a fallback when `grad != "analytic"` is requested for
// the weighted model).
static void fd_gradient(const std::vector<double>& x,
                        std::vector<double>& g, double f,
                        const XptrClosureState& s) {
  const int n = (int)x.size();
  std::vector<double> tmp = x;
  for (int i = 0; i < n; ++i) {
    const double xi = x[i];
    const double dx = std::abs(xi) * s.gradstep_rel + s.gradstep_abs;
    tmp[i] = xi + dx;
    const double f_plus = eval_value(tmp, s);
    if (s.grad_mode == GRAD_CENTRAL) {
      tmp[i] = xi - dx;
      const double f_minus = eval_value(tmp, s);
      g[i] = (f_plus - f_minus) / (2.0 * dx);
    } else {
      g[i] = (f_plus - f) / dx;
    }
    tmp[i] = xi;
    if (!std::isfinite(g[i])) g[i] = 0.0;
  }
}

} // namespace

// ===========================================================================
// create_niche_obj_ptr (refactored)
// ===========================================================================

//[[Rcpp::export]]
SEXP create_niche_obj_ptr(
    NumericMatrix env_occ,
    Nullable<NumericMatrix> env_m         = R_NilValue,
    double                  eta           = 1.0,
    std::string             likelihood    = "weighted",
    Nullable<IntegerVector> den_idx       = R_NilValue,
    Nullable<IntegerVector> kde_idx       = R_NilValue,
    Nullable<NumericVector> precomp_w_occ = R_NilValue,
    Nullable<NumericVector> precomp_w_den = R_NilValue,
    std::string             grad          = "central",
    NumericVector           gradstep      = NumericVector::create(1e-6, 1e-8)) {

  if (gradstep.size() != 2) Rcpp::stop("gradstep must have length 2");
  if (!(gradstep[0] >= 0) || !(gradstep[1] >= 0))
    Rcpp::stop("gradstep components must be non-negative.");
  if (gradstep[0] == 0.0 && gradstep[1] == 0.0)
    Rcpp::stop("gradstep cannot be both zero");

  LikType lt;
  if (likelihood == "unweighted")          lt = LIK_UNWEIGHTED;
  else if (likelihood == "weighted")       lt = LIK_WEIGHTED;
  else if (likelihood == "presence_only")  lt = LIK_PRESENCE_ONLY;
  else Rcpp::stop("Unknown likelihood type: %s", likelihood.c_str());

  GradMode gm;
  if (grad == "analytic") {
    if (lt != LIK_WEIGHTED) {
      // analytic gradient currently only implemented for the weighted model;
      // silently downgrade to central FD for the others (cheap, no regress).
      gm = GRAD_CENTRAL;
    } else {
      gm = GRAD_ANALYTIC;
    }
  } else if (grad == "central") {
    gm = GRAD_CENTRAL;
  } else if (grad == "forward") {
    gm = GRAD_FORWARD;
  } else {
    Rcpp::stop("grad must be one of 'analytic', 'central', 'forward'");
  }

  if (lt != LIK_PRESENCE_ONLY && env_m.isNull())
    Rcpp::stop("env_m must be provided for likelihood = '%s'.",
               likelihood.c_str());

  // Allocate state (shared_ptr so that the lambda capture extends lifetime
  // until the XPtr is GC'd by R).
  auto state = std::make_shared<XptrClosureState>();
  state->p           = env_occ.ncol();
  state->eta         = eta;
  state->lik_type    = lt;
  state->grad_mode   = gm;
  state->gradstep_rel = gradstep[0];
  state->gradstep_abs = gradstep[1];

  // Snapshot env_occ
  state->env_occ = Eigen::MatrixXd(env_occ.nrow(), env_occ.ncol());
  std::memcpy(state->env_occ.data(), env_occ.begin(),
              sizeof(double) * env_occ.size());

  // Snapshot env_m / build M_den + w_occ + w_den as needed
  if (lt == LIK_UNWEIGHTED) {
    NumericMatrix M(env_m);
    if (M.ncol() != state->p) Rcpp::stop("env_m must have same columns as env_occ");
    state->env_m_full = Eigen::MatrixXd(M.nrow(), M.ncol());
    std::memcpy(state->env_m_full.data(), M.begin(),
                sizeof(double) * M.size());
  } else if (lt == LIK_WEIGHTED) {
    NumericMatrix M(env_m);
    if (M.ncol() != state->p) Rcpp::stop("env_m must have same columns as env_occ");
    Eigen::Map<Eigen::MatrixXd> M_eig(
        Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(M));
    const int n_m = M.nrow();

    // Resolve denominator subset (M_den)
    if (den_idx.isNotNull()) {
      IntegerVector dv(den_idx);
      const int n_den = dv.size();
      state->M_den.resize(n_den, state->p);
      for (int i = 0; i < n_den; ++i) {
        const int row = dv[i] - 1;
        if (row < 0 || row >= n_m) Rcpp::stop("den_idx out of range");
        state->M_den.row(i) = M_eig.row(row);
      }
    } else {
      state->M_den = M_eig;
    }
    const int n_den = state->M_den.rows();

    // Resolve KDE reference subset (M_kde) - used only to compute w_occ /
    // w_den HERE if they are not precomputed by the caller.
    Eigen::MatrixXd M_kde;
    if (kde_idx.isNotNull()) {
      IntegerVector kv(kde_idx);
      const int n_kde = kv.size();
      M_kde.resize(n_kde, state->p);
      for (int i = 0; i < n_kde; ++i) {
        const int row = kv[i] - 1;
        if (row < 0 || row >= n_m) Rcpp::stop("kde_idx out of range");
        M_kde.row(i) = M_eig.row(row);
      }
    } else {
      M_kde = M_eig;
    }

    // Presence-side weights (w_occ): precompute or compute now
    if (precomp_w_occ.isNotNull()) {
      NumericVector pw(precomp_w_occ);
      if (pw.size() != env_occ.nrow())
        Rcpp::stop("precomp_w_occ length (%d) must match nrow(env_occ) (%d).",
                   (int)pw.size(), (int)env_occ.nrow());
      state->w_occ = Eigen::VectorXd(pw.size());
      for (int i = 0; i < pw.size(); ++i) state->w_occ(i) = pw[i];
    } else {
      if (state->p == 2) state->w_occ = nicher::kde_2d(state->env_occ, M_kde);
      else               state->w_occ = nicher::kde_eigen(state->env_occ, M_kde);
    }

    // Denominator-side weights (w_den)
    if (precomp_w_den.isNotNull()) {
      NumericVector pw(precomp_w_den);
      if (pw.size() != n_den)
        Rcpp::stop("precomp_w_den length (%d) must match nrow(M_den) (%d).",
                   (int)pw.size(), n_den);
      state->w_den = Eigen::VectorXd(pw.size());
      for (int i = 0; i < pw.size(); ++i) state->w_den(i) = pw[i];
    } else {
      if (state->p == 2) state->w_den = nicher::kde_2d(state->M_den, M_kde);
      else               state->w_den = nicher::kde_eigen(state->M_den, M_kde);
    }
  }

  // -------------------------------------------------------------------------
  // Build the closure
  // -------------------------------------------------------------------------
  NicheObjFunType* obj_fun = new NicheObjFunType(
    [state](const std::vector<double>& x,
            std::vector<double>& g,
            double& f) {
      const int n = (int)x.size();
      g.resize(n);
      try {
        if (state->grad_mode == GRAD_ANALYTIC && state->lik_type == LIK_WEIGHTED) {
          f = nicher::loglik_niche_math_weighted_grad_eigen(
                x.data(), n,
                state->env_occ, state->M_den, state->w_occ, state->w_den,
                state->eta,
                state->gradstep_rel, state->gradstep_abs,
                g.data());
        } else {
          f = eval_value(x, *state);
          fd_gradient(x, g, f, *state);
        }
        // Sanitize: any NaN/Inf -> abandon-direction sentinel
        if (!std::isfinite(f)) {
          f = std::numeric_limits<double>::infinity();
          std::fill(g.begin(), g.end(), 0.0);
        }
      } catch (const std::exception& e) {
        f = std::numeric_limits<double>::infinity();
        std::fill(g.begin(), g.end(), 0.0);
        REprintf("nicher XPtr ObjFun error: %s\n", e.what());
      } catch (...) {
        f = std::numeric_limits<double>::infinity();
        std::fill(g.begin(), g.end(), 0.0);
        REprintf("nicher XPtr ObjFun error: unknown C++ exception\n");
      }
    });

  return Rcpp::XPtr<NicheObjFunType>(obj_fun, true);
}
