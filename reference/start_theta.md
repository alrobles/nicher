# Starting values for niche model on math scale

Starting values for niche model on math scale

## Usage

``` r
start_theta(env_occ, skew = FALSE, skew_t = FALSE)
```

## Arguments

- env_occ:

  Data frame with environmental values at presence points.

- skew:

  Logical. If `TRUE`, append `p` zeros for the skew parameters
  (`alpha_1, ..., alpha_p`). The Gaussian-only start (`skew = FALSE`,
  the default) is used for `"presence_only"`, `"weighted"`, and
  `"ip_weighted"`; the skew start is used for `"skew_normal"` and
  `"skew_normal_weighted"`.

- skew_t:

  Logical. If `TRUE`, append `log(10)` for the skew-t degrees-of-freedom
  parameter `log_r` after the alpha block (used for `"skew_t"` and
  `"skew_t_weighted"`). Implies `skew = TRUE`.

## Value

Numeric vector of starting values for \`theta\`.

## Examples

``` r
start_theta(example_env_occ_2d)
#>      bio1WH     bio12WH      bio1WH     bio12WH             
#>   21.397746 1872.493151    1.038323    6.517823    0.000000 
start_theta(example_env_occ_2d, skew = TRUE)
#>      bio1WH     bio12WH      bio1WH     bio12WH                         
#>   21.397746 1872.493151    1.038323    6.517823    0.000000    0.000000 
#>             
#>    0.000000 
start_theta(example_env_occ_2d, skew_t = TRUE)
#>      bio1WH     bio12WH      bio1WH     bio12WH                         
#>   21.397746 1872.493151    1.038323    6.517823    0.000000    0.000000 
#>                         
#>    0.000000    2.302585 
```
