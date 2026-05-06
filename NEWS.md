# nicher 1.0.0

## CRAN release

This is the first CRAN-facing release of `nicher`.

## Model families

* `optimize_niche()` supports seven likelihood families:
  `presence_only`, `kde_bias_corrected`, `weighted`, `skew_normal`,
  `skew_normal_weighted`, `skew_t`, and `skew_t_weighted`.
* Weighted families use M-background environmental data through `env_m`.
* Non-symmetric families add skew-normal and skew-t likelihoods for
  asymmetric or heavy-tailed environmental responses.

## Model comparison and validation

* Added S3 methods for `logLik()`, `AIC()`, `BIC()`, and `nobs()`.
* Added `compare_nicher()` for side-by-side model comparison.
* Added `cv_nicher()` for k-fold and leave-one-out cross-validation of
  fitted `nicher` objects.

## Documentation

* README examples now cover all seven likelihood families.
* Added vignettes for the getting-started workflow, cross-validation,
  and non-symmetric skew-normal/skew-t models.
