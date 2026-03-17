
<!-- README.md is generated from README.Rmd. Please edit that file -->

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

``` r
library(terra)
library(nicher)

# Load raster data from extdata
stack_path <- system.file("extdata", "stack_1_12_crop.rds", package = "nicher")
example_rasters <- terra::unwrap(readRDS(stack_path))

# Extract occurrence points
data(example_occ_df)
env_occ <- terra::extract(example_rasters, example_occ_df[, c("long", "lat")])
env_occ <- env_occ[complete.cases(env_occ), -1]

# Randomly sample background points from rasters
set.seed(123) # For reproducibility
env_m <- terra::spatSample(example_rasters, size = 1000, method = "regular", na.rm = TRUE)

# Test the presence-only algorithm using optimize_niche
result <- optimize_niche(
  env_occ = env_occ,
  env_m = env_m,
  num_starts = 5,
  start_method = "uniform",
  likelihood = "presence_only"
)

# View result
print(result$best)
```
