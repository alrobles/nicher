# Model 7: Skew-t Weighted Niche

## 1 Ecological Motivation

This is the most general model in the `nicher` hierarchy. It combines
all three extensions:

1.  **Background weighting** (from Models 2/5): corrects for non-uniform
    environmental availability in \\M\\ via a KDE
    \\\hat{g}(\mathbf{x})\\.
2.  **Skewness** (from Models 4/5): allows asymmetric environmental
    tolerances via the skewness vector \\\alpha\\.
3.  **Heavy tails** (from Model 6): allows persistence at extreme
    environmental values via the degrees-of-freedom parameter \\r\\.

Use this model when both asymmetric tolerances and heavy tails are
ecologically plausible, and the observation process is confounded by
environmental availability.

## 2 Data-Generating Process

The DGP is the use-availability model with a skew-t niche kernel:

\\ p(\mathbf{x}\_i \mid \text{presence in } M) =
\frac{f_T(\mathbf{x}\_i;\mu,\Sigma,\alpha,r)\\\hat{g}(\mathbf{x}\_i)}
{\int_M
f_T(\mathbf{x};\mu,\Sigma,\alpha,r)\\\hat{g}(\mathbf{x})\\d\mathbf{x}},
\tag{1} \\

where \\f_T\\ is the skew-t density (Model 6, Eq. 2).

## 3 Mathematical Specification

### 3.1 Negative Log-Likelihood

\\ -\ell^\*(\theta) = \frac{p+r}{2}\sum\_{i=1}^{n}\log A_i -
\sum\_{i=1}^{n}\log I_i - \sum\_{i=1}^{n}\log\hat{g}(\mathbf{x}\_i) +
n\cdot\operatorname{logsumexp}\_j\\\left\[ \log\hat{g}(\mathbf{y}\_j) -
\tfrac{p+r}{2}\log A_j + \log I_j \right\] + \text{penalty}, \tag{2} \\

where:

- \\A_i = 1 + q_i/r\\ with \\q_i = \\L^{-1}(\mathbf{x}\_i - \mu)\\^2\\
  (occurrence Mahalanobis distance),
- \\A_j\\ and \\q_j\\ are the same at background point \\j\\,
- \\I_i\\ and \\I_j\\ are the Gauss–Laguerre quadrature integrals (Model
  6, Eq. 6),
- \\\hat{g}(\cdot)\\ is the precomputed KDE of the background.

### 3.2 Structure of the Denominator

The denominator integrand for background point \\j\\ is:

\\ a_j = \log\hat{g}(\mathbf{y}\_j) - \frac{p+r}{2}\log A_j + \log I_j.
\tag{3} \\

This is proportional to \\\hat{g}(\mathbf{y}\_j) \cdot
f_T(\mathbf{y}\_j)\\ evaluated in log space — the discrete approximation
to the denominator integral in Eq. (1).

### 3.3 Dropped Constants

Compared with Model 6 (Eq. 7), the terms \\\frac{n}{2}\log\|\Sigma\| +
\frac{np}{2}\log r + n\log\Gamma(r/2)\\ cancel between numerator and
denominator (they appear in both \\f_T(\mathbf{x}\_i)\\ and
\\f_T(\mathbf{y}\_j)\\) and are therefore absent from Eq. (2).

### 3.4 Relationship to All Other Models

The skew-t weighted model nests all simpler models:

| Restriction | Resulting model |
|----|----|
| \\\alpha = 0,\\r \to \infty,\\\text{no background}\\ | Model 1: Presence-only |
| \\\alpha = 0,\\r \to \infty\\ | Model 2: Weighted |
| \\\alpha = 0,\\r \to \infty,\\w = 1/\hat{g}\\ | Model 3: IP-weighted |
| \\r \to \infty,\\\text{no background}\\ | Model 4: Skew-normal |
| \\r \to \infty\\ | Model 5: Skew-normal weighted |
| \\\text{no background}\\ | Model 6: Skew-t |
| (none) | **Model 7: Skew-t weighted** |

## 4 Implementation Call Chain

    optimize_niche(likelihood = "skew_t_weighted")
      └─ niche_obj.cpp::create_niche_obj_ptr()
           └─ loglik_niche_math_skew_t_weighted_eigen()  [C++, lines 829–910]
                ├─ build_L_cov()                          [theta → mu, sigma, L_cov]
                ├─ alpha = Map(theta + 2*p + n_v, p)      [line 854]
                ├─ log_r = theta[3*p + n_v]; r = exp()    [lines 855–856]
                │
                │  // Occurrence side
                ├─ Mahalanobis q_occ, z_occ                [lines 864–868]
                ├─ For each i:
                │    A_i = 1 + q_occ(i)/r                  [line 873]
                │    sum_log_A_occ += log(A_i)             [line 874]
                │    sum_log_I_occ += log_skew_t_integral  [line 875]
                │
                │  // Denominator side
                ├─ Mahalanobis q_den, z_den                [lines 880–884]
                ├─ log_w_den = clamp(w_den).log()          [line 885]
                ├─ For each j:
                │    A_j = 1 + q_den(j)/r                  [line 889]
                │    log_I_j = log_skew_t_integral()       [line 890]
                │    log_terms(j) = log_w_den - half_kpr*log(A_j) + log_I_j  [line 891]
                ├─ logsumexp(log_terms)                    [lines 893–895]
                │
                ├─ log_w_occ = clamp(w_occ).log()          [line 897]
                ├─ prior_penalty_value(has_alpha=true)      [lines 899–900]
                └─ neg_log = half_kpr*sum_log_A_occ
                             - sum_log_I_occ
                             - log_w_occ.sum()
                             + n*lse + penalty              [lines 902–906]

### 4.1 Key Code (loglik_niche_math_cpp.cpp, lines 829–910)

The C++ kernel:

1.  Reconstructs \\L\_{\text{cov}}\\, maps \\\alpha\\ and \\\log r\\
    (lines 850–857).
2.  **Occurrence side** (lines 863–876): computes Mahalanobis distances,
    skew scores, and accumulates \\\sum_i\log A_i\\ and \\\sum_i\log
    I_i\\ via `log_skew_t_integral()`.
3.  **Denominator side** (lines 878–895): for each background point,
    computes \\A_j\\, \\\log I_j\\, and the full integrand
    \\\log\hat{g}\_j - \frac{p+r}{2}\log A_j + \log I_j\\. Applies
    logsumexp across all \\K\\ background points.
4.  Combines all terms (lines 902–906) and returns the guarded negative
    log-likelihood.

### 4.2 Computational Cost

Each evaluation requires \\2 \times Q = 64\\ calls to
`log_skew_t_integral()` per data point (once for occurrence, once for
background). For \\n = 200\\ presence points and \\K = 10{,}000\\
background points, this is \\(200 + 10{,}000) \times 32 = 326{,}400\\
\\\Phi\\ evaluations per log-likelihood call.

### 4.3 R Entry Point

``` r

optimize_niche(env_occ, env_m, likelihood = "skew_t_weighted",
               prior_log_sigma_lambda = 1.0,
               prior_alpha_lambda = 1.0)
```

## 5 When To Use

- When the **full model** is warranted: species with asymmetric,
  heavy-tailed tolerances in a region with heterogeneous environmental
  availability.
- As the **most complex candidate** in an AIC/BIC model-selection
  exercise. If the data do not justify the extra parameters (\\\alpha\\,
  \\\log r\\), simpler models will be preferred.
- When **robustness to outliers** is desired: the heavy tails of the
  skew-t down-weight extreme presence records that would strongly
  influence a Gaussian or skew-normal fit.

## 6 Limitations

1.  **Highest parameter count**. \\3p + p(p-1)/2 + 1\\ parameters. For
    \\p = 5\\, this is 26 parameters.

2.  **Slowest model**. The double loop over occurrence + background with
    \\Q = 32\\ quadrature nodes makes this the most computationally
    expensive model. Consider using `num_starts = 1` for initial
    exploration and increasing for the final fit.

3.  **Patil & Ord drift** (shared with Models 2/5). The ridge prior on
    \\\log\sigma\\ is essential.

4.  **Quadrature accuracy**. Same caveat as Model 6: \\Q = 32\\ is
    validated for \\p = 2\\.

5.  **Not all extensions are independent**. Skewness and heavy tails can
    interact in unexpected ways. Always compare with simpler nested
    models.

## References
