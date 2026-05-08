# Model 3: Inverse-Probability-Weighted Normal Niche (Horvitz--Thompson)

## 1 Ecological Motivation

When the environmental conditions available in the accessible area $`M`$
are unevenly distributed, a naive maximum-likelihood estimate of the
niche centroid $`\mu`$ is biased toward whatever environments happen to
be common in $`M`$. For example, if 70 % of $`M`$ has temperature near
15 °C, the presence-only MLE will be pulled toward 15 °C even if the
species tolerates a wider or differently centred range.

The **inverse-probability-weighted** (IPW) model removes this
environmental-availability confound by reweighting each background point
by the *inverse* of its local environmental density
$`w(\mathbf{x}) = 1/\hat{g}(\mathbf{x})`$. This converts the discrete
sum over the heterogeneous background grid into a numerical
approximation of a Lebesgue-measure integral, effectively projecting the
niche onto a uniform environmental space.

This is a Horvitz–Thompson (HT) estimator (Horvitz and Thompson (1952))
applied to the restricted likelihood framework of Jimenez et al. (2019).

## 2 Data-Generating Process

### 2.1 The Jimenez DGP

The IPW model assumes the **Jimenez et al. (2019) DGP**: presences are
drawn from the fundamental niche density
$`f(\mathbf{x};\mu,\Sigma)`$*restricted* to the set of environments that
exist in $`M`$:

``` math
p(\mathbf{x}_i \mid \text{presence in } M) =
\frac{f(\mathbf{x}_i;\mu,\Sigma)}
     {\int_M f(\mathbf{x};\mu,\Sigma)\,d\mathbf{x}}.
\tag{1}
```

Under this DGP the observation probability at environment $`\mathbf{x}`$
depends *only* on the species’ intrinsic tolerance, not on how common
$`\mathbf{x}`$ is in $`M`$. This is appropriate when:

- Sampling effort is approximately uniform across $`M`$ (e.g. systematic
  surveys).
- The analyst explicitly wants to estimate the *fundamental* niche
  (Hutchinson 1957), not the realised niche shaped by $`M`$.

### 2.2 Contrast with the Use-Availability DGP

The `"weighted"` model (Model 2) assumes instead that presence intensity
$`\propto f \cdot g`$. The two DGPs yield different estimands:

| DGP | Weight | Estimand |
|----|----|----|
| Jimenez (this model) | $`w = 1/\hat{g}`$ | Fundamental niche on uniform $`E`$-space |
| Use-availability (Model 2) | $`w = \hat{g}`$ | Niche under $`M`$-biased observation |

Neither is universally “correct”; the choice depends on the ecological
question and the data-collection protocol.

## 3 Mathematical Specification

### 3.1 Restricted Log-Likelihood

Starting from Equation 8 of Jimenez and Soberón (2022) with
$`w = 1/\hat{g}`$:

``` math
\ell(\mu,\Sigma) =
\sum_{i=1}^{n}\log f(\mathbf{x}_i)
- \sum_{i=1}^{n}\log\hat{g}(\mathbf{x}_i)
- n\,\log\!\sum_{j=1}^{K}\frac{f(\mathbf{y}_j)}{\hat{g}(\mathbf{y}_j)}.
\tag{2}
```

### 3.2 Negative Log-Likelihood (Minimised Objective)

Substituting $`f = \mathcal{N}(\mu,\Sigma)`$ and dropping the
$`(2\pi)^{-p/2}`$ constant (which cancels between numerator and
denominator):

``` math
-\ell^*(\theta) =
\frac{1}{2}\sum_{i=1}^{n} q_i
+ \sum_{i=1}^{n}\log\hat{g}(\mathbf{x}_i)
+ n\cdot\operatorname{logsumexp}_j\!\left[
  -\tfrac{1}{2}q_j - \log\hat{g}(\mathbf{y}_j)
\right],
\tag{3}
```

where $`q_i = \|L^{-1}(\mathbf{x}_i - \mu)\|^2`$ (squared Mahalanobis
distance), $`q_j`$ is the same quantity at background point
$`\mathbf{y}_j`$, and $`\hat{g}(\cdot)`$ is the Gaussian KDE evaluated
at the point.

### 3.3 Sign Convention (vs. Weighted)

| Term | `"weighted"` (Eq. 8, $`w = \hat{g}`$) | `"ip_weighted"` (Eq. 8, $`w = 1/\hat{g}`$) |
|----|----|----|
| Numerator KDE | $`-\sum_i \log\hat{g}(\mathbf{x}_i)`$ | $`+\sum_i \log\hat{g}(\mathbf{x}_i)`$ |
| Denominator integrand | $`+\log\hat{g}(\mathbf{y}_j) - \frac{1}{2}q_j`$ | $`-\frac{1}{2}q_j - \log\hat{g}(\mathbf{y}_j)`$ |

Both are instances of Equation 8; the sign inversion is **not a bug**
but the algebraic consequence of using the reciprocal weight.

### 3.4 Importance-Sampling Interpretation

The denominator
$`\sum_j \exp(-\frac{1}{2}q_j - \log\hat{g}(\mathbf{y}_j))`$ is
equivalent to the Horvitz–Thompson estimator:

``` math
\widehat{\int_M f\,d\mathbf{x}} =
\frac{1}{K}\sum_{j=1}^{K}\frac{f(\mathbf{y}_j)}{\hat{g}(\mathbf{y}_j)},
\tag{4}
```

which converges to $`\int_M f(\mathbf{x})\,d\mathbf{x}`$ as
$`K \to \infty`$ whenever $`\hat{g}(\mathbf{y}_j) > 0`$ for all $`j`$ in
the support of $`f`$ (Horvitz and Thompson (1952); Hahn (1998)).

## 4 Implementation Call Chain

    optimize_niche(likelihood = "ip_weighted")
      └─ niche_obj.cpp::create_niche_obj_ptr()
           └─ loglik_niche_math_ip_weighted_eigen()   [C++, lines 101–154]
                ├─ build_L_cov()                       [theta → mu, sigma, L_cov]
                ├─ Mahalanobis q1 (occurrences), q2 (background)
                ├─ a_j = -0.5 * q2 - log_w_den         [line 143]
                ├─ logsumexp(a)                          [lines 145–147]
                ├─ neg_log = 0.5*sum_q1 + log_w_occ + n*lse
                └─ OPTIM_PENALTY guard                   [line 152]

### 4.1 Key Code (loglik_niche_math_cpp.cpp, lines 101–154)

The C++ kernel:

1.  Reconstructs $`L_{\text{cov}}`$ from $`\theta`$ via `build_L_cov()`
    (line 123).
2.  Computes Mahalanobis distances for both occurrence and background
    points (lines 129–137).
3.  Clamps KDE weights: `w.max(MIN_KDE_WEIGHT).log()` with
    `MIN_KDE_WEIGHT = 1e-300` to prevent $`\log(0)`$ (lines 141–142).
4.  Computes the integrand
    $`a_j = -0.5 q_j - \log\hat{g}(\mathbf{y}_j)`$ (line 143).
5.  Evaluates $`\operatorname{logsumexp}(a)`$ via the max-shift trick
    (lines 145–147).
6.  Returns `0.5 * sum_q1 + log_w_occ.sum() + n * log_sum_exp` (lines
    149–151).

### 4.2 Analytic Gradient (lines 186–304)

The hybrid analytic/finite-difference gradient follows the same
structure as the weighted gradient (see Model 2):

- **$`\nabla_\mu`$**: Closed-form via
  $`\Sigma^{-1}[\text{weighted mean shift}]`$, with softmax weights
  $`\pi_j = \exp(a_j)/\sum_k\exp(a_k)`$ on the background side (lines
  263–271).
- **$`\nabla_{\log\sigma_k}`$**: Closed-form via $`U \odot V`$ products
  (lines 273–278).
- **$`\nabla_{v_k}`$**: Central finite differences on the C-vine block
  (lines 280–298).

### 4.3 R Entry Point

``` r

optimize_niche(env_occ, env_m, likelihood = "ip_weighted")
```

No ridge prior is applied by default because the IPW denominator
provides natural self-regularisation (see Section 6).

## 5 KDE Computation

The KDE $`\hat{g}`$ is computed identically to Model 2 (Scott’s rule
bandwidth, Gaussian kernel, precomputed before optimisation). The only
difference is that $`\hat{g}`$ enters the likelihood with **inverted
sign**:

- In `"weighted"`: the denominator uses $`+\log\hat{g}`$ (amplifying
  common environments).
- In `"ip_weighted"`: the denominator uses $`-\log\hat{g}`$
  (down-weighting common environments).

## 6 Self-Regularisation

Unlike the weighted model, the IPW model does **not** require a ridge
prior on $`\log\sigma`$ to prevent $`\sigma \to \infty`$. The mechanism:

1.  As $`\sigma \to \infty`$, $`f(\mathbf{y}_j) \to \text{const}`$ for
    all $`j`$.
2.  The denominator becomes
    $`\sum_j \text{const}/\hat{g}(\mathbf{y}_j)`$, which is large
    because $`1/\hat{g}`$ amplifies low-density tail points.
3.  Meanwhile the numerator $`\prod_i f(\mathbf{x}_i)`$ also flattens,
    but the denominator grows faster due to the tail amplification.
4.  The net effect is that $`-\ell^*`$*increases* as
    $`\sigma \to \infty`$, providing a natural penalty.

This self-regularisation is a property of the Horvitz–Thompson
estimator, not an ad hoc fix.

## 7 Limitations

1.  **Variance amplification**. When $`\hat{g}(\mathbf{y}_j) \approx 0`$
    (environmental deserts), the weight $`1/\hat{g}`$ can be extremely
    large, inflating the variance of the estimator. Diagnostic: check
    the effective sample size
    $`\text{ESS} = (\sum_j w_j)^2 / \sum_j w_j^2`$.

2.  **KDE bandwidth sensitivity**. Scott’s rule may under-smooth in low
    dimensions or over-smooth in high dimensions. Users should inspect
    the KDE fit.

3.  **DGP mismatch**. If presences are genuinely generated by the
    use-availability process ($`\lambda \propto f \cdot g`$), then the
    IPW estimator targets the wrong quantity. Use `"weighted"` (Model 2)
    instead.

4.  **Normality assumption**. Same as Model 1. For asymmetric niches,
    consider combining IPW with a skew-normal or skew-t kernel (not yet
    implemented in `nicher`).

## References

Hahn, J. 1998. “On the Role of the Propensity Score in Efficient
Semiparametric Estimation of Average Treatment Effects.” *Econometrica*
66 (2): 315–31. <https://doi.org/10.2307/2998560>.

Horvitz, D. G., and D. J. Thompson. 1952. “A Generalization of Sampling
Without Replacement from a Finite Universe.” *Journal of the American
Statistical Association* 47 (260): 663–85.
<https://doi.org/10.1080/01621459.1952.10483446>.

Jimenez, L., and J. Soberón. 2022. “Weighted Likelihood Estimation of
the Fundamental Niche.” *Ecological Modelling* 470: 110009.
<https://doi.org/10.1016/j.ecolmodel.2022.110009>.

Jimenez, L., J. Soberón, J. A. Christen, and D. Soto. 2019. “On the
Problem of Modeling a Fundamental Niche from Occurrence Data.”
*Ecological Modelling* 397: 74–83.
<https://doi.org/10.1016/j.ecolmodel.2019.01.020>.
