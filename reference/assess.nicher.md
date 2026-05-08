# Assess acceptance criteria for a nicher result

Evaluates the quality of a multi-start optimization by comparing
converged solutions. Returns one of four diagnostic flags:

- `"accepted_global"`:

  Multiple starts agree closely on the same high log-likelihood value —
  the global optimum is likely found.

- `"accepted_noise"`:

  Converged solutions are close but show minor numerical noise — likely
  a good solution.

- `"suggest_average"`:

  Converged solutions spread across distinct values — consider averaging
  `theta` or increasing `num_starts`.

- `"needs_more_starts"`:

  Too few converged solutions to draw conclusions — increase
  `num_starts`.

## Usage

``` r
# S3 method for class 'nicher'
assess(x, tol_gap = 0.01, tol_dist = 0.05, min_converged = 2L, ...)
```

## Arguments

- x:

  A `"nicher"` object returned by
  [`optimize_niche`](https://alrobles.github.io/nicher/reference/optimize_niche.md).

- tol_gap:

  Numeric. Relative tolerance used to decide between `"accepted_global"`
  and `"accepted_noise"`. Default 0.01.

- tol_dist:

  Numeric. Relative tolerance used to decide between `"accepted_noise"`
  and `"suggest_average"`. Default 0.05.

- min_converged:

  Integer. Minimum number of converged solutions required to avoid
  `"needs_more_starts"`. Default 2.

- ...:

  Ignored.

## Value

A named list with:

- `flag`:

  Character scalar, one of the four flags above.

- `recommendation`:

  Human-readable action recommendation.

- `gap`:

  Absolute log-likelihood gap between the best and second-best converged
  solutions.

- `rel_gap`:

  Relative gap (`gap / |best_loglik|`).

- `n_converged`:

  Number of converged solutions.

- `best_loglik`:

  Best log-likelihood value.

## Examples

``` r
# \donttest{
result <- optimize_niche(
  env_occ    = example_env_occ_2d,
  env_m      = example_env_m_2d,
  num_starts = 5L,
  likelihood = "ip_weighted"
)
diag <- assess(result)
cat("Flag          :", diag$flag, "\n")
#> Flag          : accepted_global 
cat("Recommendation:", diag$recommendation, "\n")
#> Recommendation: Global optimum likely found: multiple converged starts agree closely. 
cat("Best log-lik  :", round(diag$best_loglik, 4L), "\n")
#> Best log-lik  : -631.9638 
# }
```
