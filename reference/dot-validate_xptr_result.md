# Validate ucminfcpp::ucminf_xptr result against ucminf::ucminf.

Runs a short
[`ucminf::ucminf`](https://rdrr.io/pkg/ucminf/man/ucminf.html)
optimization from the best theta and warns if the two log-likelihoods
differ by more than `1e-3`. This guards against pointer-safety issues in
the C++ backend.

## Usage

``` r
.validate_xptr_result(
  best,
  env_occ,
  env_m,
  likelihood,
  ctrl,
  weighted_inputs = NULL,
  prior_log_sigma_center = NULL,
  prior_log_sigma_lambda = 0,
  prior_mu_center = NULL,
  prior_mu_lambda = 0,
  prior_alpha_lambda = 0,
  ...
)
```
