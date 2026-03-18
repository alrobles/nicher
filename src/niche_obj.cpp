// src/niche_obj.cpp
//
// NicheObjFun: A C++ functor wrapping the nicher negative log-likelihood
// in a format compatible with ucminfcpp's ucminf::ObjFun (std::function).
//
// The external pointer returned by create_niche_obj_ptr() can be passed
// directly to ucminfcpp::ucminf_xptr() from R, eliminating R interpreter
// overhead on every function/gradient evaluation.
//
// Supports all three likelihood types implemented in nicher:
//   "unweighted"    — loglik_niche_chol_cpp
//   "weighted"      — loglik_niche_weighted_integrated_cpp
//   "presence_only" — loglik_niche_presence_only_cpp

// [[Rcpp::depends(RcppEigen)]]
#include "nicher_types.h"
#include <functional>
#include <vector>
#include <cmath>
#include <memory>
#include <string>

using namespace Rcpp;

// ---------------------------------------------------------------------------
// Forward declarations of C++ functions defined in other files of this package
// ---------------------------------------------------------------------------

NumericMatrix cvine_cholesky(NumericVector v, int d, double eta);

double loglik_niche_chol_cpp(NumericVector mu,
                              NumericMatrix L,
                              NumericMatrix env_occ,
                              NumericMatrix env_m);

double loglik_niche_presence_only_cpp(NumericVector mu,
                                       NumericMatrix L,
                                       NumericMatrix env_occ);

double loglik_niche_weighted_integrated_cpp(NumericVector mu,
                                             NumericMatrix L,
                                             NumericMatrix env_occ,
                                             NumericMatrix env_m,
                                             Nullable<IntegerVector> den_idx,
                                             Nullable<IntegerVector> kde_idx,
                                             Nullable<NumericVector> precomp_w_den,
                                             bool neg);

// ---------------------------------------------------------------------------
// Type alias compatible with ucminf::ObjFun in ucminfcpp
// (ucminf::ObjFun is itself defined as this same std::function alias)
// ---------------------------------------------------------------------------
using NicheObjFunType = std::function<void(const std::vector<double>&,
                                            std::vector<double>&,
                                            double&)>;

// ---------------------------------------------------------------------------
// NicheObjFunData: stores all model data needed by the functor
// Uses Eigen matrices for C++-heap storage (safe from R's garbage collector)
// ---------------------------------------------------------------------------
struct NicheObjFunData {
    Eigen::MatrixXd env_occ;          // presence-point environmental matrix (n_occ x p)
    Eigen::MatrixXd env_m;            // background environmental matrix    (n_m   x p)
    double          eta;              // LKJ prior shape parameter
    int             p;                // number of environmental dimensions
    int             lik_type;         // 0=unweighted, 1=weighted, 2=presence_only

    // Optional subsampling data for the weighted likelihood
    std::vector<int>    den_idx_vec;       // 1-based denominator indices (R convention)
    std::vector<int>    kde_idx_vec;       // 1-based KDE-reference indices
    std::vector<double> precomp_w_den_vec; // precomputed denominator KDE weights
    bool has_den_idx;
    bool has_kde_idx;
    bool has_precomp_w;

    // ------------------------------------------------------------------
    // Compute the negative log-likelihood at parameter vector x.
    //
    // Parameter encoding (same as in loglik_niche_math_cpp):
    //   x[0 .. p-1]      : mu  (centroid)
    //   x[p .. 2p-1]     : log_sigma  (log of standard deviations)
    //   x[2p .. end]     : v  (C-vine partial correlation parameters)
    // ------------------------------------------------------------------
    double compute(const std::vector<double>& x) const {
        int n = static_cast<int>(x.size());

        // Decode mu and sigma = exp(log_sigma)
        NumericVector mu(p), sigma(p);
        for (int i = 0; i < p; ++i) {
            mu[i]    = x[i];
            sigma[i] = std::exp(x[p + i]);
        }

        // Decode C-vine parameters v
        int nv = n - 2 * p;
        NumericVector v(nv);
        for (int i = 0; i < nv; ++i) v[i] = x[2 * p + i];

        // Build Cholesky factor: L_cov = diag(sigma) * L_corr
        NumericMatrix L_corr = cvine_cholesky(v, p, eta);
        NumericMatrix L_cov(p, p);
        for (int j = 0; j < p; ++j)
            for (int k = 0; k <= j; ++k)
                L_cov(j, k) = sigma[j] * L_corr(j, k);

        // Wrap stored Eigen data as Rcpp matrices (zero-copy via Map)
        NumericMatrix occ_m(env_occ.rows(), env_occ.cols());
        Eigen::Map<Eigen::MatrixXd>(occ_m.begin(),
                                    env_occ.rows(),
                                    env_occ.cols()) = env_occ;

        if (lik_type == 2) {
            // Presence-only: no background data needed
            return loglik_niche_presence_only_cpp(mu, L_cov, occ_m);
        }

        NumericMatrix m_m(env_m.rows(), env_m.cols());
        Eigen::Map<Eigen::MatrixXd>(m_m.begin(),
                                    env_m.rows(),
                                    env_m.cols()) = env_m;

        if (lik_type == 0) {
            // Unweighted
            return loglik_niche_chol_cpp(mu, L_cov, occ_m, m_m);
        }

        // Weighted (lik_type == 1)
        Nullable<IntegerVector> di = R_NilValue;
        Nullable<IntegerVector> ki = R_NilValue;
        Nullable<NumericVector> pw = R_NilValue;

        if (has_den_idx) {
            IntegerVector dv(den_idx_vec.begin(), den_idx_vec.end());
            di = dv;
        }
        if (has_kde_idx) {
            IntegerVector kv(kde_idx_vec.begin(), kde_idx_vec.end());
            ki = kv;
        }
        if (has_precomp_w) {
            NumericVector pv(precomp_w_den_vec.begin(), precomp_w_den_vec.end());
            pw = pv;
        }

        return loglik_niche_weighted_integrated_cpp(mu, L_cov, occ_m, m_m,
                                                     di, ki, pw, true);
    }
};

// ---------------------------------------------------------------------------
// create_niche_obj_ptr
//
// Creates a C++ objective function (wrapped in a std::function compatible
// with ucminf::ObjFun) that evaluates the nicher negative log-likelihood
// and its central-difference gradient entirely in C++, then returns it
// as an Rcpp::XPtr for use with ucminfcpp::ucminf_xptr().
//
// @param env_occ        Numeric matrix of environmental values at presences.
// @param env_m          Numeric matrix of background environmental values.
//                       May be NULL when likelihood = "presence_only".
// @param eta            LKJ prior shape parameter (default 1.0).
// @param likelihood     One of "unweighted", "weighted", "presence_only".
// @param den_idx        Integer vector of 1-based indices for denominator
//                       subsample (weighted model only). NULL = use all.
// @param kde_idx        Integer vector of 1-based indices for KDE reference
//                       subsample (weighted model only). NULL = use all.
// @param precomp_w_den  Precomputed denominator KDE weights (optional).
// @param gradstep       Length-2 numeric vector: relative and absolute step
//                       sizes for central finite differences.
//                       Default c(1e-6, 1e-8) matches the ucminf default.
//                       Step for parameter i: |x_i| * gradstep[0] + gradstep[1]
//
// @return An external pointer (externalptr) wrapping a heap-allocated
//         NicheObjFunType (std::function).  Pass to ucminfcpp::ucminf_xptr().
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
SEXP create_niche_obj_ptr(NumericMatrix env_occ,
                           Nullable<NumericMatrix> env_m       = R_NilValue,
                           double                  eta         = 1.0,
                           std::string             likelihood  = "unweighted",
                           Nullable<IntegerVector> den_idx     = R_NilValue,
                           Nullable<IntegerVector> kde_idx     = R_NilValue,
                           Nullable<NumericVector> precomp_w_den = R_NilValue,
                           NumericVector           gradstep    = NumericVector::create(1e-6, 1e-8)) {

    if (gradstep.size() != 2)
        Rcpp::stop("gradstep must have length 2");

    // Validate likelihood type
    int lt;
    if      (likelihood == "unweighted")    lt = 0;
    else if (likelihood == "weighted")      lt = 1;
    else if (likelihood == "presence_only") lt = 2;
    else
        Rcpp::stop("Unknown likelihood type '%s'. "
                   "Use 'unweighted', 'weighted', or 'presence_only'.",
                   likelihood.c_str());

    if (lt != 2 && env_m.isNull())
        Rcpp::stop("env_m must be provided for likelihood '%s'.",
                   likelihood.c_str());

    // Build data object on the heap (owned by the shared_ptr)
    auto data = std::make_shared<NicheObjFunData>();
    data->eta      = eta;
    data->p        = env_occ.ncol();
    data->lik_type = lt;

    // Deep-copy matrices into Eigen storage (safe from R's GC)
    data->env_occ = Eigen::Map<Eigen::MatrixXd>(
        env_occ.begin(), env_occ.nrow(), env_occ.ncol());

    if (lt != 2) {
        NumericMatrix m(env_m);
        data->env_m = Eigen::Map<Eigen::MatrixXd>(
            m.begin(), m.nrow(), m.ncol());
    }

    // Store optional subsampling data
    data->has_den_idx  = den_idx.isNotNull();
    data->has_kde_idx  = kde_idx.isNotNull();
    data->has_precomp_w = precomp_w_den.isNotNull();

    if (data->has_den_idx) {
        IntegerVector dv(den_idx);
        data->den_idx_vec.assign(dv.begin(), dv.end());
    }
    if (data->has_kde_idx) {
        IntegerVector kv(kde_idx);
        data->kde_idx_vec.assign(kv.begin(), kv.end());
    }
    if (data->has_precomp_w) {
        NumericVector pv(precomp_w_den);
        data->precomp_w_den_vec.assign(pv.begin(), pv.end());
    }

    // Gradient finite-difference step sizes.
    // Step for parameter i: |x_i| * gradstep_rel + gradstep_abs
    // This matches the formula used by ucminf (same defaults: 1e-6, 1e-8).
    const double gs_rel = gradstep[0];
    const double gs_abs = gradstep[1];

    // Allocate the std::function (NicheObjFunType) on the heap.
    // The shared_ptr to data is captured so the data lives as long as the
    // function object does.
    NicheObjFunType* obj_fun = new NicheObjFunType(
        [data, gs_rel, gs_abs]
        (const std::vector<double>& x,
         std::vector<double>&       g,
         double&                    f)
        {
            // Evaluate objective
            f = data->compute(x);

            // Central finite differences for gradient
            int n = static_cast<int>(x.size());
            g.resize(n);
            for (int i = 0; i < n; ++i) {
                double xi = x[i];
                double dx = std::abs(xi) * gs_rel + gs_abs;
                std::vector<double> xp = x, xm = x;
                xp[i] = xi + dx;
                xm[i] = xi - dx;
                g[i]  = (data->compute(xp) - data->compute(xm)) / (2.0 * dx);
            }
        }
    );

    // Return as an external pointer with automatic deletion via finalizer
    return Rcpp::XPtr<NicheObjFunType>(obj_fun, true);
}
