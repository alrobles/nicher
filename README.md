
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nicher

<!-- badges: start -->

<!-- badges: end -->

The goal of nicher is to create ecological niche models based on an
ellipse model.

## Installation

You can install the development version of nicher from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("alrobles/nicher")
```

## Example

``` r
library(nicher)

# Built-in 2D hummingbird environmental datasets
env_occ <- example_env_occ_2d   # presences
env_m   <- example_env_m_2d     # background (accessibility area M)

# Run multi-start optimization using Sobol design
# breadth = 0.1 → starting quantiles c(0.1, 0.5, 0.9)
result <- optimize_niche(
  env_occ    = env_occ,
  env_m      = env_m,
  num_starts = 20L,
  breadth    = 0.1,
  likelihood = "unweighted"
)

# Print a concise summary of the optimization result
print(result)
#> -- nicher optimization result --
#>   Likelihood : unweighted
#>   Starts     : 20 (converged: 20 )
#>   Best loglik: -3.141593
#>   Convergence: 1

# Assess acceptance criteria
diag <- assess(result)
cat("Flag          :", diag$flag, "\n")
cat("Recommendation:", diag$recommendation, "\n")
cat("Best log-lik  :", round(diag$best_loglik, 4L), "\n")
cat("Gap           :", round(diag$gap, 6L), "\n")
#> Flag          : accepted_global
#> Recommendation: Global optimum likely found: multiple converged starts agree closely.
#> Best log-lik  : -3.1416
#> Gap           : 0.000012

# Extract the best parameter vector
theta_best <- result$best$theta
```

The four possible `assess()` flags are:

| Flag | Meaning |
|------|---------|
| `accepted_global` | Multiple starts agree → global optimum likely found |
| `accepted_noise` | Good solution, minor numerical noise across starts |
| `suggest_average` | Wide spread → consider averaging `theta` or more starts |
| `needs_more_starts` | Too few converged solutions → increase `num_starts` |

