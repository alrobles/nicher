# Non-Symmetric Models: Skew-Normal, Skew-t, and NCST Families

## Motivation

Symmetric Gaussian models assume that a species responds equally to
deviations above and below its environmental optimum. Many species
violate this assumption: heat tolerance may extend further than cold
tolerance, or precipitation responses may be right-skewed.

`nicher` provides three non-symmetric likelihood families that relax the
symmetry assumption while retaining the ellipsoidal geometry:

| Family | Parameters | Tail weight | Background correction |
|----|----|----|----|
| `skew_normal` | $`3p + p(p-1)/2`$ | Gaussian | No |
| `skew_t` | $`3p + p(p-1)/2 + 1`$ | Heavy | No |
| `ncst` | $`3p + p(p-1)/2 + 1`$ | Heavy | No |
| `skew_normal_weighted` | $`3p + p(p-1)/2`$ | Gaussian | Yes (KDE) |
| `skew_t_weighted` | $`3p + p(p-1)/2 + 1`$ | Heavy | Yes (KDE) |
| `ncst_weighted` | $`3p + p(p-1)/2 + 1`$ | Heavy | Yes (KDE) |

## 1 The skew-normal family

The skew-normal distribution (Azzalini and Capitanio 1999) adds a vector
$`\alpha \in \mathbb{R}^p`$ of skewness parameters to the multivariate
normal. When $`\alpha = 0`$, the model reduces to the symmetric
Gaussian.

### Fitting

``` r

library(nicher)
set.seed(42)

fit_sn <- optimize_niche(
  env_occ    = example_env_occ_2d,
  env_m      = NULL,
  likelihood = "skew_normal",
  num_starts = 5L,
  breadth    = 0.1,
  control    = list(maxeval = 1000L),
  verbose    = FALSE
)
fit_sn
#> -- nicher optimization result --
#>   Likelihood : skew_normal 
#>   Starts     : 5 (converged: 5 )
#>   Best loglik: -668.1408 
#>   Convergence: 1 
#>   eta        : 1 
#>   Variables  : bio1WH, bio12WH
```

The estimated skewness parameters indicate the direction and magnitude
of asymmetry:

``` r

theta <- fit_sn$best$theta
p <- ncol(example_env_occ_2d)
n_v <- p * (p - 1L) / 2L
alpha <- theta[(2L * p + n_v + 1L):(3L * p + n_v)]
names(alpha) <- colnames(example_env_occ_2d)
alpha
#>      bio1WH     bio12WH 
#> -0.02343198  3.27622204
```

## 2 The skew-t family

The skew-t model (Branco and Dey 2001) extends the skew-normal by
replacing Gaussian tails with Student-$`t`$ tails controlled by a
degrees-of-freedom parameter $`r > 0`$:

``` math
T = \mu + X\sqrt{r/Y}, \qquad X \sim \text{SN}_p(0, \Sigma, \alpha),
\quad Y \sim \chi^2_r.
```

As $`r \to \infty`$, the skew-t converges to the skew-normal.

### Fitting

``` r

set.seed(42)

fit_st <- optimize_niche(
  env_occ    = example_env_occ_2d,
  env_m      = NULL,
  likelihood = "skew_t",
  num_starts = 5L,
  breadth    = 0.1,
  control    = list(maxeval = 1000L),
  verbose    = FALSE
)
fit_st
#> -- nicher optimization result --
#>   Likelihood : skew_t 
#>   Starts     : 5 (converged: 5 )
#>   Best loglik: -718.787 
#>   Convergence: 1 
#>   eta        : 1 
#>   Variables  : bio1WH, bio12WH
```

The estimated degrees of freedom:

``` r

log_r <- fit_st$best$theta[length(fit_st$best$theta)]
cat("Estimated r =", round(exp(log_r), 2), "\n")
#> Estimated r = 170.02
```

## 3 The NCST family

The Non-Central Skew-$`t`$ (NCST) places the location parameter
$`\xi`$*inside* the skew-normal component before chi-squared scaling:

``` math
T = \frac{X}{\sqrt{Y/r}}, \qquad X \sim \text{SN}_p(\xi, \Omega, \alpha),
\quad Y \sim \chi^2_r.
```

This allows degrees of freedom to interact with the location-scale
structure, providing extra flexibility compared to the standard skew-t.

### Fitting

``` r

set.seed(42)

fit_ncst <- optimize_niche(
  env_occ    = example_env_occ_2d,
  env_m      = NULL,
  likelihood = "ncst",
  num_starts = 5L,
  breadth    = 0.1,
  control    = list(maxeval = 1000L),
  verbose    = FALSE
)
fit_ncst
#> -- nicher optimization result --
#>   Likelihood : ncst 
#>   Starts     : 5 (converged: 5 )
#>   Best loglik: -719.3935 
#>   Convergence: 1 
#>   eta        : 1 
#>   Variables  : bio1WH, bio12WH
```

## 4 Model comparison

Compare all non-symmetric families using AIC and BIC:

``` r

cmp <- compare_nicher(
  skew_normal = fit_sn,
  skew_t      = fit_st,
  ncst        = fit_ncst
)
print(cmp[, c("model", "likelihood", "AIC", "dAIC", "BIC", "dBIC")])
#>         model  likelihood      AIC     dAIC      BIC     dBIC
#> 1 skew_normal skew_normal 1348.135   0.0000 1364.168   0.0000
#> 2      skew_t      skew_t 1451.464 103.3291 1469.788 105.6195
#> 3        ncst        ncst 1453.603 105.4687 1471.927 107.7592
```

## 5 Weighted variants

Each non-symmetric family has a weighted counterpart that corrects for
sampling bias using a KDE of the background environment. Supply `env_m`
and use the `_weighted` suffix:

``` r

set.seed(42)

fit_snw <- optimize_niche(
  env_occ    = example_env_occ_2d,
  env_m      = example_env_m_2d,
  likelihood = "skew_normal_weighted",
  num_starts = 5L,
  breadth    = 0.1,
  control    = list(maxeval = 500L),
  verbose    = FALSE
)
fit_snw
#> -- nicher optimization result --
#>   Likelihood : skew_normal_weighted 
#>   Starts     : 6 (converged: 6 )
#>   Best loglik: -652.4987 
#>   Convergence: 1 
#>   eta        : 1 
#>   Variables  : bio1WH, bio12WH
```

## 6 When to use non-symmetric models

- **Skew-normal**: When the species has clearly asymmetric environmental
  tolerances but Gaussian-like tails.
- **Skew-t**: When the species tolerates extreme environmental values
  more than a Gaussian would predict (heavy tails).
- **NCST**: When additional flexibility in how tail heaviness interacts
  with location is needed. Compare with skew-t via AIC/BIC.
- **Weighted variants**: When sampling bias in the background
  environment must be corrected simultaneously with asymmetry.

Use
[`compare_nicher()`](https://alrobles.github.io/nicher/reference/compare_nicher.md)
or
[`cv_nicher()`](https://alrobles.github.io/nicher/reference/cv_nicher.md)
to select among families for a given dataset.

## References

Azzalini, A., and A. Capitanio. 1999. “Statistical Applications of the
Multivariate Skew Normal Distribution.” *Journal of the Royal
Statistical Society: Series B* 61 (3): 579–602.
<https://doi.org/10.1111/1467-9868.00194>.

Branco, M. D., and D. K. Dey. 2001. “A General Class of Multivariate
Skew-Elliptical Distributions.” *Journal of Multivariate Analysis* 79
(1): 99–113. <https://doi.org/10.1006/jmva.2000.1960>.
