# Why the Weighted Model Struggles on the Vicugna Dataset

## Introduction

The `nicher` package ships with a real-world dataset for *Vicugna
vicugna* (the wild ancestor of the alpaca), containing 545 occurrence
records and 10 000 background points from the species’ accessible area
(M) in the central Andes, characterised by two bioclimatic variables:
annual mean temperature (bio1, in degrees Celsius) and annual
precipitation (bio12, in mm).

This vignette explains **why the weighted likelihood family produces
worse results than the presence-only family on this dataset**, and what
statistical properties of the data drive the difference. Understanding
this case study is important for practitioners choosing among the seven
likelihood families available in `nicher`.

``` r

library(nicher)
data(example_vicugna)

env_occ <- example_vicugna$env_occ
env_m   <- example_vicugna$env_m
```

## The core issue: multimodal M with a concentrated niche

The vicugna occupies a narrow band in environmental space: most
occurrences cluster tightly around bio1 \\\approx\\ 5–6 degrees Celsius
and bio12 \\\approx\\ 100 mm. The accessible area (M), however, spans a
much wider and **multimodal** range of climates.

``` r

cat("Occurrences (n =", nrow(env_occ), ")\n")
#> Occurrences (n = 545 )
print(summary(env_occ))
#>       bio1             bio12     
#>  Min.   :-0.4523   Min.   :   8  
#>  1st Qu.: 3.9153   1st Qu.:  68  
#>  Median : 5.6048   Median : 111  
#>  Mean   : 6.1184   Mean   : 215  
#>  3rd Qu.: 8.0380   3rd Qu.: 274  
#>  Max.   :15.9072   Max.   :1061

cat("\nBackground M (n =", nrow(env_m), ")\n")
#> 
#> Background M (n = 10000 )
print(summary(env_m))
#>       bio1            bio12       
#>  Min.   :-7.926   Min.   :   4.0  
#>  1st Qu.: 4.462   1st Qu.: 105.0  
#>  Median : 6.916   Median : 236.0  
#>  Mean   : 6.876   Mean   : 387.2  
#>  3rd Qu.: 8.916   3rd Qu.: 723.0  
#>  Max.   :19.403   Max.   :1455.0
```

The key diagnostic is the distribution of bio12 (precipitation) in M:

``` r

d_m   <- density(env_m$bio12)
d_occ <- density(env_occ$bio12)

# Count density peaks
peaks_m <- which(diff(sign(diff(d_m$y))) == -2) + 1L
cat("Density peaks in M (bio12):", length(peaks_m), "\n")
#> Density peaks in M (bio12): 3
cat("Peak locations:", round(d_m$x[peaks_m], 1), "mm\n")
#> Peak locations: 100.8 853.2 1327.7 mm
```

The M distribution has multiple modes: a primary mode around 100 mm
(arid Altiplano), a secondary mode near 850 mm (eastern Andean slopes),
and a minor shoulder at higher precipitation. This multimodality has
direct consequences for the weighted likelihood.

## The presence-only model: a clean baseline

The presence-only family estimates \\(\mu, \Sigma)\\ by maximising the
log-likelihood: \\\ell\_{\mathrm{PO}}(\theta) = \sum\_{i=1}^{n} \log
f(x_i \mid \theta)\\ where \\f\\ is the multivariate Gaussian density.
Since it only uses the occurrence data, the PO estimator reduces to the
sample mean and covariance, which are well-defined for the unimodal,
concentrated vicugna niche.

``` r

set.seed(42)
fit_po <- optimize_niche(
  env_occ    = env_occ,
  env_m      = NULL,
  num_starts = 20L,
  likelihood = "presence_only",
  seed       = 42L
)
cat("PO loglik:", fit_po$best$loglik, "\n")
#> PO loglik: -4116.357
cat("PO convergence:", fit_po$best$convergence, "\n")
#> PO convergence: 1
cat("PO AIC:", AIC(fit_po), "\n")
#> PO AIC: 8242.715

p <- ncol(env_occ)
cat("PO mu:", round(fit_po$best$theta[seq_len(p)], 2), "\n")
#> PO mu: 6.12 214.96
```

The PO estimates are stable and close to the empirical mean: bio1
\\\approx\\ 6.1, bio12 \\\approx\\ 215.

## The weighted model: denominator instability

The weighted likelihood (Jimenez et al. 2022) introduces a correction
for the background: \\\ell\_{\mathrm{W}}(\theta) = \sum\_{i=1}^{n} \log
\frac{f(x_i \mid \theta)}{\hat g(x_i)} - \log \sum\_{j=1}^{m}
\frac{f(x_j^{(M)} \mid \theta)}{\hat g(x_j^{(M)})}\\ where \\\hat g\\ is
a kernel density estimate of the background distribution. The
denominator (the integral approximation) sums over background points,
weighted inversely by the KDE.

``` r

set.seed(42)
fit_w <- optimize_niche(
  env_occ    = env_occ,
  env_m      = env_m,
  num_starts = 20L,
  likelihood = "weighted",
  seed       = 42L
)
cat("W loglik:", fit_w$best$loglik, "\n")
#> W loglik: -5062.32
cat("W convergence:", fit_w$best$convergence, "\n")
#> W convergence: 1
cat("W AIC:", AIC(fit_w), "\n")
#> W AIC: 10120.82

cat("W mu:", round(fit_w$best$theta[seq_len(p)], 2), "\n")
#> W mu: -0.62 235.01
```

The weighted model’s \\\mu\\ estimate for bio1 is pulled far from the
occurrence centre, and the log-sigma values are inflated. **The model
fits worse by every metric**: higher (worse) AIC, and a \\\mu\\ that
does not correspond to the observed niche.

### Why this happens

Three interacting factors create the problem:

1.  **Multimodal KDE denominator.** When \\\hat g(x)\\ has multiple
    modes, the inverse weights \\1/\hat g(x)\\ amplify points in the
    low-density valleys between modes. Background points near bio12
    \\\approx\\ 400–600 mm (the trough between the 100 mm and 850 mm
    peaks) receive disproportionately high weight, pulling the
    denominator integral and distorting the gradient landscape.

2.  **Scale mismatch between niche and M.** The occurrences occupy a
    small, concentrated region of E-space, but M extends over a far
    wider volume. The denominator integral involves evaluating the
    Gaussian \\f(x \mid \theta)\\ over a 10 000-point grid that mostly
    lies far from \\\mu\\. Small changes in \\\Sigma\\ cause large
    changes in the denominator sum, creating a nearly flat or multimodal
    objective.

3.  **Heavy right tail in bio12.** The precipitation variable has a
    mean/median ratio of \\\approx 1.9\\ in the occurrences and
    \\\approx 1.6\\ in M, indicating strong right-skewness. A symmetric
    Gaussian kernel in the weighted likelihood cannot represent this
    shape faithfully. The optimizer compensates by inflating
    \\\sigma\_{\mathrm{bio12}}\\ and shifting \\\mu\\.

## Model comparison

``` r

cmp <- compare_nicher(
  presence_only = fit_po,
  weighted      = fit_w
)
#> Warning: Mixing fits with and without `env_m` (likelihood families
#> presence_only, weighted). AIC/BIC values are still computed but interpret with
#> care: weighted families normalise the likelihood by an env_m integral, so their
#> loglik is on a different scale than presence_only / skew_normal / skew_t.
print(cmp[, c("model", "likelihood", "AIC", "dAIC", "BIC", "dBIC")])
#>           model    likelihood       AIC    dAIC       BIC    dBIC
#> 1 presence_only presence_only  8242.715    0.00  8264.219    0.00
#> 2      weighted      weighted 10120.825 1878.11 10142.329 1878.11
```

The presence-only model is decisively preferred by both AIC and BIC.

## When does the weighted model help?

The weighted correction is valuable when:

- **M is unimodal** (or approximately so) — the KDE weights are stable.
- **The niche is broad relative to M** — the denominator integral is
  well-conditioned.
- **Sampling bias is strong** — occurrences are geographically clustered
  in a way that does not reflect environmental preference. The weighted
  likelihood corrects for this by down-weighting over-sampled
  environments.

For vicugna, sampling is relatively uniform across the species’ range
and the niche is narrow, so the correction introduces more noise than
signal.

## Recommendations for practitioners

1.  **Always fit presence-only first** as a baseline. The `warm_start`
    mechanism in
    [`optimize_niche()`](https://alrobles.github.io/nicher/reference/optimize_niche.md)
    does this automatically for weighted fits.

2.  **Compare models with AIC/BIC** using
    [`compare_nicher()`](https://alrobles.github.io/nicher/reference/compare_nicher.md).
    A weighted model with worse AIC than presence-only indicates that
    the background correction is not helping.

3.  **Consider skew-normal families** when bio12 or other variables show
    strong asymmetry. The `skew_normal` family adds shape parameters
    \\\alpha\\ that can capture the right tail without inflating
    \\\Sigma\\.

4.  **Inspect the M distribution** before choosing a weighted family. If
    `density(env_m$variable)` shows multiple peaks, the KDE-based
    correction may be unstable.

5.  **Use
    [`cv_nicher()`](https://alrobles.github.io/nicher/reference/cv_nicher.md)**
    to tune regularisation strength when the weighted model is needed.
    Cross-validation can select penalty parameters (`prior_mu_lambda`,
    `prior_log_sigma_lambda`) that stabilise the weighted fit.

## References

Jimenez, L., Soberon, J., Christen, J. A., & Soto, D. (2022). Dealing
with critical uncertainty in ecological niche models. *Ecological
Modelling*, 468, 109823.
<https://doi.org/10.1016/j.ecolmodel.2021.109823>
