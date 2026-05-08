# Model 5: Skew-Normal Weighted Niche

## 1 Ecological Motivation

This model combines two extensions:

1.  **Skewness** (from Model 4): the niche can be asymmetric along each
    environmental axis.
2.  **Background weighting** (from Model 2): the non-uniform
    distribution of environments in \\M\\ is modelled via a KDE
    \\\hat{g}(\mathbf{x})\\.

The result is a model that estimates an **asymmetric fundamental niche**
while accounting for the fact that presence records are more likely in
environments that are common in \\M\\ (the use-availability DGP).

## 2 Data-Generating Process

The DGP is the same as Model 2, but with a skew-normal kernel replacing
the Gaussian:

\\ p(\mathbf{x}\_i \mid \text{presence in } M) =
\frac{f\_{\text{SN}}(\mathbf{x}\_i;\mu,\Sigma,\alpha)\\\hat{g}(\mathbf{x}\_i)}
{\int_M
f\_{\text{SN}}(\mathbf{x};\mu,\Sigma,\alpha)\\\hat{g}(\mathbf{x})\\d\mathbf{x}},
\tag{1} \\

where \\f\_{\text{SN}}\\ is the skew-normal density (Model 4, Eq. 1).

## 3 Mathematical Specification

### 3.1 Negative Log-Likelihood

\\ -\ell^\*(\theta) = \frac{1}{2}\sum\_{i=1}^{n} q_i -
\sum\_{i=1}^{n}\log\Phi(z_i) -
\sum\_{i=1}^{n}\log\hat{g}(\mathbf{x}\_i) +
n\cdot\operatorname{logsumexp}\_j\\\left\[ \log\hat{g}(\mathbf{y}\_j) +
\log\Phi(z_j) - \tfrac{1}{2}q_j \right\] + \text{penalty}, \tag{2} \\

where:

- \\q_i = \\L^{-1}(\mathbf{x}\_i - \mu)\\^2\\ (squared Mahalanobis
  distance at presence point \\i\\),
- \\q_j\\ is the same at background point \\j\\,
- \\z_i = \alpha^\top\operatorname{diag}(\sigma)^{-1}(\mathbf{x}\_i -
  \mu)\\ is the skew score at presence point \\i\\,
- \\z_j\\ is the skew score at background point \\j\\,
- \\\hat{g}(\cdot)\\ is the Gaussian KDE of the background,
- the penalty term is as in Model 4, Eq. (4).

### 3.2 Structure of the Denominator

The denominator integrand (inside the logsumexp) has three log-scale
terms:

\\ a_j =
\underbrace{\log\hat{g}(\mathbf{y}\_j)}\_{\text{availability}} +
\underbrace{\log\Phi(z_j)}\_{\text{skewness}} -
\underbrace{\tfrac{1}{2}q_j}\_{\text{Mahalanobis}}. \tag{3} \\

This is exactly the log of \\\hat{g}(\mathbf{y}\_j) \cdot \Phi(z_j)
\cdot \exp(-\frac{1}{2}q_j)\\, which is proportional to \\\hat{g} \cdot
f\_{\text{SN}}\\ evaluated at background point \\j\\ — the discrete
approximation to the denominator integral in Eq. (1).

### 3.3 Relationship to Models 2 and 4

| Component       | Model 2 (weighted) | Model 4 (skew-normal) | Model 5 (this) |
|-----------------|--------------------|-----------------------|----------------|
| Mahalanobis     | \\\checkmark\\     | \\\checkmark\\        | \\\checkmark\\ |
| \\\log\Phi(z)\\ | —                  | \\\checkmark\\        | \\\checkmark\\ |
| KDE weighting   | \\\checkmark\\     | —                     | \\\checkmark\\ |
| Background data | Required           | Not used              | Required       |

## 4 Implementation Call Chain

    optimize_niche(likelihood = "skew_normal_weighted")
      └─ niche_obj.cpp::create_niche_obj_ptr()
           └─ loglik_niche_math_skew_normal_weighted_eigen()  [C++, lines 590–656]
                ├─ build_L_cov()                               [theta → mu, sigma, L_cov]
                ├─ alpha = Map(theta + 2*p + n_v, p)           [line 615]
                ├─ Mahalanobis q1 (occ), q2 (background)
                ├─ sum_log_phi_occ = sum_log_pnorm(occ)        [line 624]
                ├─ z_den = skew_z(diff_den, sigma, alpha)      [line 632]
                ├─ log_phi_den = log_pnorm(z_den)              [lines 633–634]
                ├─ a_j = log_w_den + log_phi_den - 0.5*q2      [line 640]
                ├─ logsumexp(a)                                  [lines 642–644]
                ├─ prior_penalty_value(has_alpha=true)           [lines 646–647]
                └─ neg_log = 0.5*sum_q1 - sum_log_phi_occ
                             - log_w_occ + n*lse + penalty       [lines 649–653]

### 4.1 Key Code (loglik_niche_math_cpp.cpp, lines 590–656)

The C++ kernel:

1.  Reconstructs \\L\_{\text{cov}}\\ and maps \\\alpha\\ (lines
    611–615).
2.  Computes occurrence-side Mahalanobis distances and
    \\\sum_i\log\Phi(z_i)\\ (lines 620–624).
3.  Computes background-side Mahalanobis distances \\q_j\\ and skew
    scores \\z_j\\ (lines 626–634).
4.  Evaluates the denominator integrand \\a_j = \log\hat{g}\_j +
    \log\Phi(z_j) - \frac{1}{2}q_j\\ (line 640).
5.  Applies the logsumexp trick (lines 642–644).
6.  Adds the ridge + alpha penalty (lines 646–647).
7.  Returns the guarded negative log-likelihood (lines 649–655).

### 4.2 R Entry Point

``` r

optimize_niche(env_occ, env_m, likelihood = "skew_normal_weighted",
               prior_log_sigma_lambda = 1.0,
               prior_alpha_lambda = 1.0)
```

## 5 When To Use

- When both **asymmetric tolerances** and **non-uniform environmental
  availability** are expected.
- As a model-selection alternative to `"weighted"` (Model 2): fit both
  and compare via AIC/BIC. If the skewness parameters are close to zero,
  the simpler Model 2 is preferred.
- For species with known thermal or precipitation asymmetries in regions
  where the accessible area has heterogeneous environments.

## 6 Limitations

1.  **Parameter count**. Adds \\p\\ skewness parameters over Model 2.
    Requires larger sample sizes for reliable estimation.

2.  **Patil & Ord drift** (shared with Model 2). The ridge prior on
    \\\log\sigma\\ is essential; do not set \\\lambda = 0\\.

3.  **Unimodal**. The skew-normal is still unimodal. For bimodal niches,
    mixture extensions are needed (not implemented).

4.  **Numerical stability**. The product \\\Phi(z_j) \cdot \hat{g}\_j\\
    can be very small for background points far from the niche centre.
    The logsumexp trick (lines 642–644) handles this by working entirely
    in log space.

## References
