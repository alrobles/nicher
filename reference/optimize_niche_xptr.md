# Optimize niche model using a compiled XPtr backend

Runs ucminfcpp::ucminf_xptr() over a compiled C++ objective function
created by create_niche_obj_ptr(). Supports single-start and multi-start
optimization. Ensures that all starting vectors are numeric doubles and
finite.

## Usage

``` r
optimize_niche_xptr(
  start = NULL,
  xptr,
  control = ucminfcpp::ucminf_control(grad = "central", gradstep = c(1e-06, 1e-08),
    maxeval = 200),
  multi_start = FALSE
)
```

## Arguments

- start:

  Numeric vector (single start) or list of numeric vectors
  (multi-start).

- xptr:

  External pointer created by create_niche_obj_ptr().

- control:

  List of control parameters for ucminfcpp::ucminf_control().

- multi_start:

  Logical; TRUE if multiple starts are supplied.

## Value

A list with:

- par — best parameter vector

- value — negative log-likelihood

- conv — convergence code

- all_results — (multi-start only) table of all runs

## Examples

``` r
# \donttest{
## optimize_niche() uses optimize_niche_xptr() internally;
## call it directly for full control over each optimization run:
fit <- optimize_niche(
  env_occ    = example_env_occ_3d,
  env_m      = example_env_m_3d,
  num_starts = 3L,
  likelihood = "ip_weighted"
)
fit$best$value
#> NULL
# }
```
