# Package index

## All functions

- [`AIC(`*`<nicher>`*`)`](https://alrobles.github.io/nicher/reference/AIC.nicher.md)
  : Akaike information criterion for a nicher fit

- [`BIC(`*`<nicher>`*`)`](https://alrobles.github.io/nicher/reference/BIC.nicher.md)
  : Bayesian information criterion for a nicher fit

- [`assess()`](https://alrobles.github.io/nicher/reference/assess.md) :
  Assess acceptance criteria for an optimization result

- [`assess(`*`<nicher>`*`)`](https://alrobles.github.io/nicher/reference/assess.nicher.md)
  : Assess acceptance criteria for a nicher result

- [`compare_nicher()`](https://alrobles.github.io/nicher/reference/compare_nicher.md)
  : Compare nicher fits side-by-side via AIC and BIC

- [`cv_nicher()`](https://alrobles.github.io/nicher/reference/cv_nicher.md)
  : k-fold (or leave-one-out) cross-validation for a nicher fit

- [`cvine_cholesky()`](https://alrobles.github.io/nicher/reference/cvine_cholesky.md)
  : Build Cholesky factor of a correlation matrix from C‑vine partial
  correlations

- [`example_env_m_2d`](https://alrobles.github.io/nicher/reference/example_env_m_2d.md)
  : Samples of environmental data from M hypothesis to estimate negative
  log likelihood from Abeillia abeillei presence points. This is a
  hummingbird example. A dataset with two variables containing points
  contains 2 bioclimatic variables

- [`example_env_m_3d`](https://alrobles.github.io/nicher/reference/example_env_m_3d.md)
  : Samples points from M hypothesis to estimate negative log likelihood
  from Abeillia abeillei presence points. This is a hummingbird example.
  A dataset with three columns containing extracted information from 3
  bioclimatic variables

- [`example_env_occ_2d`](https://alrobles.github.io/nicher/reference/example_env_occ_2d.md)
  : Species occurrence points to estimate negative log likelihood from
  Abeillia abeillei presence points. This is a hummingbird example. A
  dataset with two variables containing points contains 2 bioclimatic
  variables

- [`example_env_occ_3d`](https://alrobles.github.io/nicher/reference/example_env_occ_3d.md)
  : Species occurrence points to estimate negative log likelihood from
  Abeillia abeillei presence points. This is a hummingbird example. A
  dataset with three variables containing points contains 3 bioclimatic
  variables

- [`example_mu_vec`](https://alrobles.github.io/nicher/reference/example_mu_vec.md)
  : Example of a vector of the center of an ellipsoid from two
  environmental variables. vector of length 2 corresponding to the
  centroid of an ellipsoid

- [`example_occ_df`](https://alrobles.github.io/nicher/reference/example_occ_df.md)
  : Species occurrence points from Abeillia abeillei presence points
  after download and clean from GBIF. This is a hummingbird example. A
  dataset with three variables. Contains scientific name, longitude and
  latitude.

- [`example_s_mat`](https://alrobles.github.io/nicher/reference/example_s_mat.md)
  : Example of a covariance matrix of an ellipsoid from two
  environmental variables. The 2 x 2 matrix corresponded to a positive
  semi-definite matrix. In two dimensions encodes the rotation
  (orientation) and scaling of an ellipse.

- [`example_vicugna`](https://alrobles.github.io/nicher/reference/example_vicugna.md)
  : Samples of environmental data from M hypothesis to estimate negative
  log likelihood from Vicugna vicugna.

- [`geom_nicher_background()`](https://alrobles.github.io/nicher/reference/geom_nicher_background.md)
  : Background environment layer

- [`geom_nicher_ellipse()`](https://alrobles.github.io/nicher/reference/geom_nicher_ellipse.md)
  :

  Niche-ellipse layer derived from a fitted `nicher` model

- [`geom_nicher_isosuitability()`](https://alrobles.github.io/nicher/reference/geom_nicher_isosuitability.md)
  :

  Iso-suitability contour layer derived from a fitted `nicher` model

- [`geom_nicher_occ()`](https://alrobles.github.io/nicher/reference/geom_nicher_occ.md)
  : Occurrence-points layer

- [`get_ellipsoid_pars()`](https://alrobles.github.io/nicher/reference/get_ellipsoid_pars.md)
  : Get ellipsoid parameters. A function to compute average and the
  inverse of covariance matrix from environmental data

- [`habitat_suitability()`](https://alrobles.github.io/nicher/reference/habitat_suitability.md)
  :

  Tiled habitat-suitability map from an environmental terra stack

- [`kde_gaussian()`](https://alrobles.github.io/nicher/reference/kde_gaussian.md)
  : Multivariate Gaussian KDE with fixed Scott bandwidth

- [`logLik(`*`<nicher>`*`)`](https://alrobles.github.io/nicher/reference/logLik.nicher.md)
  : Log-likelihood of a fitted nicher model

- [`loglik_niche()`](https://alrobles.github.io/nicher/reference/loglik_niche.md)
  : Negative log likelihood of an ellipsoid corrected with environmental
  combinations which come from the area of study (M)

- [`loglik_niche_math_cpp()`](https://alrobles.github.io/nicher/reference/loglik_niche_math_cpp.md)
  : Negative log-likelihood (M-restricted, math scale, Cholesky version)

- [`loglik_niche_math_ip_weighted()`](https://alrobles.github.io/nicher/reference/loglik_niche_math_ip_weighted.md)
  : Negative log-likelihood (inverse-probability-weighted normal, math
  scale)

- [`loglik_niche_math_ip_weighted_integrated()`](https://alrobles.github.io/nicher/reference/loglik_niche_math_ip_weighted_integrated.md)
  : Negative log-likelihood (inverse-probability-weighted normal,
  integrated C++)

- [`loglik_niche_math_presence_only()`](https://alrobles.github.io/nicher/reference/loglik_niche_math_presence_only.md)
  : Negative log-likelihood (presence-only, math scale)

- [`niche_ip_weighted()`](https://alrobles.github.io/nicher/reference/niche_ip_weighted.md)
  : Fit inverse-probability-weighted (IPW) normal niche model

- [`niche_presence_only()`](https://alrobles.github.io/nicher/reference/niche_presence_only.md)
  : Fit presence-only Gaussian niche model

- [`nobs(`*`<nicher>`*`)`](https://alrobles.github.io/nicher/reference/nobs.nicher.md)
  : Number of observations used to fit a nicher model

- [`optimize_niche()`](https://alrobles.github.io/nicher/reference/optimize_niche.md)
  : Optimize niche model log-likelihood with multi-start Sobol design

- [`optimize_niche_xptr()`](https://alrobles.github.io/nicher/reference/optimize_niche_xptr.md)
  : Optimize niche model using a compiled XPtr backend

- [`predict(`*`<nicher>`*`)`](https://alrobles.github.io/nicher/reference/predict.nicher.md)
  :

  Habitat-suitability raster from a fitted `nicher` object

- [`print(`*`<nicher>`*`)`](https://alrobles.github.io/nicher/reference/print.nicher.md)
  : Print a nicher object

- [`start_theta()`](https://alrobles.github.io/nicher/reference/start_theta.md)
  : Starting values for niche model on math scale

- [`start_theta_multiple()`](https://alrobles.github.io/nicher/reference/start_theta_multiple.md)
  : Generate multiple starting points for niche model optimization
