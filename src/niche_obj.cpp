// src/niche_obj.cpp — VERSION 2026-03-21
#include "nicher_types.h"
#include <functional>
#include <vector>
#include <cmath>
#include <memory>
#include <string>
using namespace Rcpp;

// forward declarations
NumericMatrix cvine_cholesky(NumericVector v, int d, double eta);

double loglik_niche_chol_cpp(
    NumericVector mu,
    NumericMatrix L,
    NumericMatrix env_occ,
    NumericMatrix env_m
);

double loglik_niche_presence_only_cpp(
    NumericVector mu,
    NumericMatrix L,
    NumericMatrix env_occ
);

double loglik_niche_weighted_integrated_cpp(
    NumericVector mu,
    NumericMatrix L,
    NumericMatrix env_occ,
    NumericMatrix env_m,
    Nullable<IntegerVector> den_idx,
    Nullable<IntegerVector> kde_idx,
    Nullable<NumericVector> precomp_w_den,
    bool neg
);


// ========================================================================
// Functor container
// ========================================================================

using NicheObjFunType =
  std::function<void(const std::vector<double>&,
                     std::vector<double>&,
                     double&)>;

struct NicheObjFunData {

  Eigen::MatrixXd env_occ;  // OWNED COPY
  Eigen::MatrixXd env_m;    // OWNED COPY (except in presence-only)
  double eta;
  int p;
  int lik_type; // 0=unweighted,1=weighted,2=presence_only

  std::vector<int> den_idx_vec;
  std::vector<int> kde_idx_vec;
  std::vector<double> precomp_w_den_vec;

  bool has_den;
  bool has_kde;
  bool has_precomp;

  // =====================================================================
  // Compute: negative log-likelihood evaluation
  // =====================================================================
  double compute(const std::vector<double>& x) const {

    // dimensional checks
    int n = x.size();
    int expected = 2*p + p*(p-1)/2;
    if (n != expected)
      stop("Parameter length mismatch in compute().");

    // unpack parameters
    NumericVector mu(p), sigma(p);
    for(int i=0;i<p;i++){
      mu[i] = x[i];
      sigma[i] = std::exp(x[p+i]);
    }
    int nv = n - 2*p;
    NumericVector v(nv);
    for(int i=0;i<nv;i++)
      v[i] = x[2*p + i];

    // build covariance Cholesky
    NumericMatrix Lcorr = cvine_cholesky(v, p, eta);
    NumericMatrix Lcov(p,p);
    for(int j=0;j<p;j++)
      for(int k=0;k<=j;k++)
        Lcov(j,k) = sigma[j] * Lcorr(j,k);

    // build OWNED local copies for loglik calls
    NumericMatrix occ_m(env_occ.rows(), env_occ.cols());
    std::memcpy(occ_m.begin(),
                env_occ.data(),
                sizeof(double)*env_occ.size());

    NumericMatrix m_m;
    if(lik_type != 2){
      m_m = NumericMatrix(env_m.rows(), env_m.cols());
      std::memcpy(m_m.begin(),
                  env_m.data(),
                  sizeof(double)*env_m.size());
    }

    // route
    if(lik_type == 2){
      return loglik_niche_presence_only_cpp(mu, Lcov, occ_m);
    }

    if(lik_type == 0){
      return loglik_niche_chol_cpp(mu, Lcov, occ_m, m_m);
    }

    // weighted
    Nullable<IntegerVector> di = R_NilValue;
    Nullable<IntegerVector> ki = R_NilValue;
    Nullable<NumericVector> pw = R_NilValue;

    if(has_den) {
      di = wrap(IntegerVector(den_idx_vec.begin(), den_idx_vec.end()));
    }
    if(has_kde){
      ki = wrap(IntegerVector(kde_idx_vec.begin(), kde_idx_vec.end()));
    }
    if(has_precomp){
      pw = wrap(NumericVector(precomp_w_den_vec.begin(),
                              precomp_w_den_vec.end()));
    }

    return loglik_niche_weighted_integrated_cpp(
      mu, Lcov, occ_m, m_m,
      di, ki, pw, true
    );
  }

};


// ========================================================================
// create_niche_obj_ptr — FINAL VERSION
// ========================================================================

//[[Rcpp::export]]
SEXP create_niche_obj_ptr(NumericMatrix env_occ,
                          Nullable<NumericMatrix> env_m = R_NilValue,
                          double eta = 1.0,
                          std::string likelihood = "unweighted",
                          Nullable<IntegerVector> den_idx = R_NilValue,
                          Nullable<IntegerVector> kde_idx = R_NilValue,
                          Nullable<NumericVector> precomp_w_den = R_NilValue,
                          NumericVector gradstep = NumericVector::create(1e-6,1e-8))
{
  if(gradstep.size() != 2)
    stop("gradstep must have length 2");

  int lt;
  if(likelihood == "unweighted") lt = 0;
  else if(likelihood == "weighted") lt = 1;
  else if(likelihood == "presence_only") lt = 2;
  else stop("Unknown likelihood type");

  if(lt != 2 && env_m.isNull())
    stop("env_m must be provided for this likelihood.");

  // allocate container
  auto data = std::make_shared<NicheObjFunData>();
  data->eta = eta;
  data->p   = env_occ.ncol();
  data->lik_type = lt;

  // --- DEEP COPY env_occ ---
  data->env_occ = Eigen::MatrixXd(env_occ.nrow(), env_occ.ncol());
  std::memcpy(
    data->env_occ.data(),
    env_occ.begin(),
    sizeof(double) * env_occ.size()
  );

  // --- DEEP COPY env_m (if exists) ---
  if(lt != 2){
    NumericMatrix M(env_m);
    if(M.ncol() != data->p)
      stop("env_m must have same columns as env_occ");

    data->env_m = Eigen::MatrixXd(M.nrow(), M.ncol());
    std::memcpy(
      data->env_m.data(),
      M.begin(),
      sizeof(double) * M.size()
    );
  }

  // copy indices
  data->has_den = den_idx.isNotNull();
  data->has_kde = kde_idx.isNotNull();
  data->has_precomp = precomp_w_den.isNotNull();

  if(data->has_den){
    IntegerVector dv(den_idx);
    data->den_idx_vec.assign(dv.begin(), dv.end());
  }

  if(data->has_kde){
    IntegerVector kv(kde_idx);
    data->kde_idx_vec.assign(kv.begin(), kv.end());
  }

  if(data->has_precomp){
    NumericVector pv(precomp_w_den);
    data->precomp_w_den_vec.assign(pv.begin(), pv.end());
  }

  double gs_rel = gradstep[0];
  double gs_abs = gradstep[1];

  if(gs_rel < 0 || gs_abs < 0)
    stop("gradstep components must be non-negative.");

  if(gs_rel==0 && gs_abs==0)
    stop("gradstep cannot be both zero");

  // ===================================================
  // FUNCTOR
  // ===================================================
  NicheObjFunType* obj_fun = new NicheObjFunType(
    [data, gs_rel, gs_abs]
    (const std::vector<double>& x,
     std::vector<double>& g,
     double& f)
    {
      f = data->compute(x);

      int n = x.size();
      g.resize(n);

      std::vector<double> tmp = x;

      for(int i=0;i<n;i++){
        double xi = x[i];
        double dx = std::abs(xi)*gs_rel + gs_abs;

        tmp[i] = xi + dx;
        double fp = data->compute(tmp);

        tmp[i] = xi - dx;
        double fm = data->compute(tmp);

        tmp[i] = xi;

        g[i] = (fp - fm) / (2.0 * dx);
      }
    }
  );

  return Rcpp::XPtr<NicheObjFunType>(obj_fun, true);
}
