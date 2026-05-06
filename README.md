
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nicher

<!-- badges: start -->

<!-- badges: end -->

`nicher` estimates ecological niche models with ellipsoidal
environmental responses under the *M* hypothesis. It includes symmetric
Gaussian models, background-weighted models, and non-symmetric
skew-normal and skew-t models for niches with asymmetric or heavy-tailed
environmental responses.

## Installation

You can install the development version of nicher from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("alrobles/nicher")
```

## Model families

The main interface is `optimize_niche()`. The current likelihood
families are:

| Family                   | Uses `env_m`? | Shape                             | When to use it                                               |
| ------------------------ | ------------: | --------------------------------- | ------------------------------------------------------------ |
| `"presence_only"`        |            No | Symmetric Gaussian                | Baseline model using occurrences only                        |
| `"kde_bias_corrected"`   |           Yes | Symmetric Gaussian                | Legacy nicher 2.x weighted behavior                          |
| `"weighted"`             |           Yes | Symmetric Gaussian                | Paper-faithful weighted likelihood with ridge regularization |
| `"skew_normal"`          |            No | Non-symmetric skew-normal         | Presence-only niches with asymmetric responses               |
| `"skew_normal_weighted"` |           Yes | Non-symmetric skew-normal         | Asymmetric niches with background correction                 |
| `"skew_t"`               |            No | Non-symmetric heavy-tailed skew-t | Asymmetric niches with heavy tails                           |
| `"skew_t_weighted"`      |           Yes | Non-symmetric heavy-tailed skew-t | Heavy-tailed asymmetric niches with background correction    |

## A small runnable example

The examples below use small subsets of the bundled hummingbird data so
that the README can be knitted quickly. For real analyses, use all
available occurrence/background data and increase `num_starts` and
`control$maxeval`.

``` r
library(nicher)

env_occ <- example_env_occ_2d[seq_len(40), ]
env_m   <- example_env_m_2d[seq_len(500), ]

doc_control <- list(maxeval = 150L)

fit_model <- function(likelihood) {
  uses_background <- likelihood %in% c(
    "kde_bias_corrected",
    "weighted",
    "skew_normal_weighted",
    "skew_t_weighted"
  )

  optimize_niche(
    env_occ = env_occ,
    env_m = if (uses_background) env_m else NULL,
    num_starts = 3L,
    likelihood = likelihood,
    seed = 1L,
    m_subsample = nrow(env_m),
    m_kde_subsample = nrow(env_m),
    control = doc_control,
    verbose = FALSE
  )
}
```

Fit the symmetric Gaussian families:

``` r
fits_symmetric <- list(
  presence_only = fit_model("presence_only"),
  kde_bias_corrected = fit_model("kde_bias_corrected"),
  weighted = fit_model("weighted")
)

vapply(fits_symmetric, function(x) x$best$loglik, numeric(1))
#>      presence_only kde_bias_corrected           weighted
#>          -344.7938          -230.3105          -249.3001
```

Fit the non-symmetric skew-normal and skew-t families:

``` r
fits_nonsymmetric <- list(
  skew_normal = fit_model("skew_normal"),
  skew_normal_weighted = fit_model("skew_normal_weighted"),
  skew_t = fit_model("skew_t"),
  skew_t_weighted = fit_model("skew_t_weighted")
)

vapply(fits_nonsymmetric, function(x) x$best$loglik, numeric(1))
#>          skew_normal skew_normal_weighted               skew_t
#>            -369.9437            -247.9566            -397.7071
#>      skew_t_weighted
#>            -240.2834
```

Each result is an S3 `nicher` object:

``` r
cat("Weighted best log-likelihood:", fits_symmetric$weighted$best$loglik, "\n")
#> Weighted best log-likelihood: -249.3001

diag <- assess(fits_symmetric$weighted)
cat("Flag          :", diag$flag, "\n")
#> Flag          : accepted_global
cat("Recommendation:", diag$recommendation, "\n")
#> Recommendation: Global optimum likely found: multiple converged starts agree closely.
```

## Parameter layout

All models use an unconstrained parameter vector `theta`. The first
blocks are always the centroid `mu`, `log(sigma)`, and C-vine
correlation parameters. Skew models append a skewness vector `alpha`;
skew-t models append one more `log_r` degrees-of-freedom parameter.

``` r
theta_summary <- data.frame(
  model = names(c(fits_symmetric, fits_nonsymmetric)),
  n_parameters = vapply(c(fits_symmetric, fits_nonsymmetric), function(x) {
    length(x$best$theta)
  }, integer(1)),
  row.names = NULL
)

theta_summary
#>                  model n_parameters
#> 1        presence_only            5
#> 2   kde_bias_corrected            5
#> 3             weighted            5
#> 4          skew_normal            7
#> 5 skew_normal_weighted            7
#> 6               skew_t            8
#> 7      skew_t_weighted            8
```

## Comparing related families

`compare_nicher()` reports log-likelihood, AIC, BIC, and Akaike weights.
Compare models fit to the same occurrence data. Weighted and
presence-only families use different normalizing assumptions, so compare
within related groups and interpret mixed comparisons carefully.

``` r
compare_nicher(
  skew_normal_weighted = fits_nonsymmetric$skew_normal_weighted,
  skew_t_weighted = fits_nonsymmetric$skew_t_weighted
)
#>                  model           likelihood    loglik df nobs      AIC     dAIC
#> 1      skew_t_weighted      skew_t_weighted -237.2819  8   40 490.5639  0.00000
#> 2 skew_normal_weighted skew_normal_weighted -246.2550  7   40 506.5099 15.94603
#>        BIC     dBIC   weight_AIC convergence
#> 1 504.0749  0.00000 0.9996554807           1
#> 2 518.3321 14.25715 0.0003445193           1
```

## Projecting suitability

Use `predict()` with a `terra` raster whose layer names match the fitted
environmental variables. The projected surface is a standardized
ellipsoidal suitability surface with values in `(0, 1]`.

``` r
env_m_matrix <- as.matrix(env_m)

template <- terra::rast(
  nrows = 40,
  ncols = 40,
  xmin = min(env_m_matrix[, 1L]),
  xmax = max(env_m_matrix[, 1L]),
  ymin = min(env_m_matrix[, 2L]),
  ymax = max(env_m_matrix[, 2L])
)

xy <- terra::crds(template, df = TRUE)
env_rast <- c(template, template)
terra::values(env_rast[[1L]]) <- xy[, 1L]
terra::values(env_rast[[2L]]) <- xy[, 2L]
names(env_rast) <- fits_symmetric$weighted$var_names

suitability <- predict(fits_symmetric$weighted, env_rast)
terra::global(suitability, c("min", "max"), na.rm = TRUE)
#>                      min       max
#> suitability 0.0005104779 0.9973846
```

## More examples

See `vignettes/nicher-intro.Rmd` for the package workflow,
`vignettes/non_symmetric_models.Rmd` for skew-normal and skew-t models,
and `vignettes/cross_validation.Rmd` for penalty tuning and model
selection.

The four possible `assess()` flags are:

| Flag                | Meaning                                                      |
| ------------------- | ------------------------------------------------------------ |
| `accepted_global`   | Multiple starts agree; a global optimum is likely            |
| `accepted_noise`    | Good solution with minor numerical noise                     |
| `suggest_average`   | Wide spread; consider averaging `theta` or increasing starts |
| `needs_more_starts` | Too few converged solutions; increase `num_starts`           |
