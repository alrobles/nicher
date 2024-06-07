
<!-- README.md is generated from README.Rmd. Please edit that file -->

![](man/figures/nicher-logo.png)

# nicher

<!-- badges: start -->
<!-- badges: end -->

The goal of nicher is to create ecological niche models based on an
ellipse model

## Installation

You can install the development version of nicher from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("alrobles/nicher")
```

## Example

First we get the shapefile of an accesibility area (M):

``` r
library(nicher)
library(terra)
#> Warning: package 'terra' was built under R version 4.2.3
#> terra 1.7.39
```

``` r


M_path <- system.file("extdata", "Mshp_test.rds", package="nicher")
Mshp <- terra::unwrap(readr::read_rds(M_path))
plot(Mshp)
```

<img src="man/figures/README-example-1.png" width="100%" />

Then we get environmental variables to model:

``` r

stack_path <- system.file("extdata", "stack_1_12_crop.rds", package="nicher")
# 2 variables
stack_1_12 <- terra::unwrap(readr::read_rds(stack_path))
plot(stack_1_12)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

``` r
# 3 variables
stack_1_12_19 <- get_example_data("stack_1_12_19")
plot(stack_1_12_19)
```

<img src="man/figures/README-unnamed-chunk-2-2.png" width="100%" />

We get the parameters of the ellipse:

``` r
pars_2var <-  get_ENM_par(rawSpOccPnts, stack_1_12, Mshp, method = "mahalanobis")

pars_3var <-  get_ENM_par(rawSpOccPnts, stack_1_12_19, Mshp, method = "mahalanobis")
```

Then we plot the ellipse (2 vars case):
<img src="man/figures/README-2varsPars-1.png" width="100%" />

Then we plot the ellipse (3 vars case):
<img src="man/figures/README-3varsPars-1.png" width="100%" />

We predict the suitability given environmental data (2 vars case):
<img src="man/figures/README-predict2vars-1.png" width="100%" />

Finally we predict the suitability given environmental data (3 vars case
All world):
<img src="man/figures/README-predict3vars-1.png" width="100%" />
