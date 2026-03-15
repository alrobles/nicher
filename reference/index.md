# Package index

## All functions

- [`create_mask_niche()`](http://alrobles.github.io/nicher/reference/create_mask_niche.md)
  : Crear una máscara para parámetros de modelo de nicho gaussiano

- [`ellipse()`](http://alrobles.github.io/nicher/reference/ellipse.md) :
  Ellipse constructor

- [`ellipse_par_example`](http://alrobles.github.io/nicher/reference/ellipse_par_example.md)
  :

  An ellipse parameter object example. Is a list with two elements. \#'
  @format A list with two elements. One vector and one matrix:

  mu

  :   mu parameters of the ellipse

  S

  :   Covariance parameters of the ellipse

- [`ellipse_points()`](http://alrobles.github.io/nicher/reference/ellipse_points.md)
  : Ellipse points

- [`get_A_matrix()`](http://alrobles.github.io/nicher/reference/get_A_matrix.md)
  : Calculates the inverse of covariance matrix from environmental
  variables

- [`get_ENM_par()`](http://alrobles.github.io/nicher/reference/get_ENM_par.md)
  : get_ENM_par

- [`get_density()`](http://alrobles.github.io/nicher/reference/get_density.md)
  : Environmental density from presence points

- [`get_ellip_par()`](http://alrobles.github.io/nicher/reference/get_ellip_par.md)
  : Get ellipsoid parameters. A function to compute average and the
  inverse of covariance matrix from environmental data

- [`get_env_var()`](http://alrobles.github.io/nicher/reference/get_env_var.md)
  : get_env_var Function to extract environmental values from presence
  points

- [`get_example_data()`](http://alrobles.github.io/nicher/reference/get_example_data.md)
  : Function to return mammal shapefiles of Rodentia order

- [`get_mahalanobis()`](http://alrobles.github.io/nicher/reference/get_mahalanobis.md)
  : Mahalanobis distance

- [`get_negative_log()`](http://alrobles.github.io/nicher/reference/get_negative_log.md)
  : Get the negative log-likelihood value from quadratic forms of
  environmental information from presence points and a sample from an M
  hypothesis

- [`get_negloglike_optim_par()`](http://alrobles.github.io/nicher/reference/get_negloglike_optim_par.md)
  : get_negloglike_optim_par Function to optime negloglike function
  given parameters from presence points and a sample of environmental
  points from a M hypothesis. Returns a list of

- [`get_negloglike_optimr_par()`](http://alrobles.github.io/nicher/reference/get_negloglike_optimr_par.md)
  : get_negloglike_optim_par Function to optime negloglike function
  given parameters from presence points and a sample of environmental
  points from a M hypothesis. Internally runs Returns a list of
  parameters

- [`get_optim_par()`](http://alrobles.github.io/nicher/reference/get_optim_par.md)
  : get_optim_par Get the parameters for the optimization of the
  Mahalanobis ellipse

- [`get_suitability()`](http://alrobles.github.io/nicher/reference/get_suitability.md)
  : get_suitability Function that calculates the log(suitability)

- [`loglik_bio_niche()`](http://alrobles.github.io/nicher/reference/loglik_bio_niche.md)
  : Log-verosimilitud para modelo de nicho gaussiano (escala biológica)

- [`loglik_math_niche()`](http://alrobles.github.io/nicher/reference/loglik_math_niche.md)
  : Log-verosimilitud para modelo de nicho gaussiano en escala
  matemática

- [`math_to_bio_niche()`](http://alrobles.github.io/nicher/reference/math_to_bio_niche.md)
  : Convertir parámetros de escala matemática a biológica para modelo de
  nicho

- [`negloglike_multivariable()`](http://alrobles.github.io/nicher/reference/negloglike_multivariable.md)
  : Negative log likelyhood

- [`plot(`*`<ellipse>`*`)`](http://alrobles.github.io/nicher/reference/plot.ellipse.md)
  : Plot method for "ellipse" class

- [`predict(`*`<ellipse>`*`)`](http://alrobles.github.io/nicher/reference/predict.ellipse.md)
  : Predict method for "ellipse" class

- [`print(`*`<ellipse>`*`)`](http://alrobles.github.io/nicher/reference/print.ellipse.md)
  : Print method for "ellipse" class

- [`quad()`](http://alrobles.github.io/nicher/reference/quad.md) : quad
  Function that calculates quadratic terms

- [`rawSpOccPnts`](http://alrobles.github.io/nicher/reference/rawSpOccPnts.md)
  :

  Raw species ocurrence points from Abeillia abeillei presence points
  after download and clean from GBIF. This is a hummingbird example. A
  dataset with three variables. Contains scientific name, longitude and
  latitude \#' @format A data frame with 73 rows and 3 variables:

  species

  :   Species name

  long

  :   Decimal longitide geographical cooordinate, in degrees

  lat

  :   Decimal latitude geographical cooordinate, in degrees

- [`samMPts`](http://alrobles.github.io/nicher/reference/samMPts.md) :

  Samples points from M hypothesis to estimate negative log likelihood
  from Abeillia abeillei presence points. This is a hummingbird example.
  A dataset with two variables containing points contains 2 bioclimatic
  variables \#' @format A data frame with 10000 rows and 2 variables:

  bio1WH

  :   temperature, in C

  bio12WH

  :   precipitaion, in mm

- [`samMPts_1_12_19`](http://alrobles.github.io/nicher/reference/samMPts_1_12_19.md)
  :

  Samples points from M hypothesis to estimate negative log likelihood
  from Abeillia abeillei presence points. This is a hummingbird example.
  A dataset with two variables containing points contains 2 bioclimatic
  variables \#' @format A data frame with 10000 rows and 2 variables:

  bio1WH

  :   Yemperature, in C

  bio12WH

  :   Precipitaion, in mm

  bio19WH

  :   Precipitation of Coldest Quarter, in mm

- [`sam_polyM()`](http://alrobles.github.io/nicher/reference/sam_polyM.md)
  : Sample environmental points from M hypothesis

- [`spOccPnts`](http://alrobles.github.io/nicher/reference/spOccPnts.md)
  :

  Species ocurrence points to estimate negative log likelihood from
  Abeillia abeillei presence points. This is a hummingbird example. A
  dataset with two variables containing points contains 2 bioclimatic
  variables \#' @format A data frame with 73 rows and 2 variables:

  bio1WH

  :   temperature, in C

  bio12WH

  :   precipitaion, in mm

- [`spOccPnts_1_12_19`](http://alrobles.github.io/nicher/reference/spOccPnts_1_12_19.md)
  :

  Species ocurrence points to estimate negative log likelihood from
  Abeillia abeillei presence points. This is a hummingbird example. A
  dataset with two variables containing points contains 2 bioclimatic
  variables \#' @format A data frame with 73 rows and 2 variables:

  bio1WH

  :   temperature, in C

  bio12WH

  :   precipitaion, in mm

  bio19WH

  :   Precipitation of Coldest Quarter, in mm
