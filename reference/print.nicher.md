# Print a nicher object

Displays a concise summary of a `"nicher"` optimization result,
including the likelihood type, number of starts, convergence count, and
the best log-likelihood achieved.

## Usage

``` r
# S3 method for class 'nicher'
print(x, ...)
```

## Arguments

- x:

  A `"nicher"` object returned by
  [`optimize_niche`](https://alrobles.github.io/nicher/reference/optimize_niche.md).

- ...:

  Ignored.

## Value

Invisibly returns `x`.

## Examples

``` r
# \donttest{
result <- optimize_niche(
  env_occ    = example_env_occ_2d,
  env_m      = example_env_m_2d,
  num_starts = 5L,
  likelihood = "ip_weighted"
)
print(result)
#> -- nicher optimization result --
#>   Likelihood : ip_weighted 
#>   Starts     : 6 (converged: 6 )
#>   Best loglik: -631.9638 
#>   Convergence: 1 
#>   eta        : 1 
#>   Variables  : bio1WH, bio12WH 
# }
```
