# Model 6: Skew-t Niche (Branco & Dey 2001)

## 1 Ecological Motivation

Some species have environmental tolerances that are not only asymmetric
but also *heavy-tailed*: they can persist at extreme environmental
values with higher probability than a Gaussian or skew-normal model
would predict. This is common for generalist species or species at range
edges where occasional extreme events (heat waves, droughts) are
tolerated.

The **skew-t** model extends the skew-normal (Model 4) by adding a
degrees-of-freedom parameter \\r \> 0\\ that controls tail heaviness:

- As \\r \to \infty\\, the skew-t converges to the skew-normal.
- For small \\r\\ (e.g. \\r = 3\\–\\5\\), the tails are substantially
  heavier than Gaussian.

The distribution follows the noncentral skew-t (NCST) of Branco and Dey
(2001).

## 2 Mathematical Specification

### 2.1 Density

The skew-t density is derived as a scale mixture of skew-normals:

\\ T = \mu + X\sqrt{r/Y}, \qquad X \sim \text{SN}\_p(0, \Sigma, \alpha),
\quad Y \sim \chi^2_r. \tag{1} \\

The resulting density at \\\mathbf{t}\\ is:

\\ f_T(\mathbf{t}) =
C(r,\Sigma)\\\left(\frac{2}{A(\mathbf{t})}\right)^{(p+r)/2}\\I(\mathbf{t}),
\tag{2} \\

where:

\\ A(\mathbf{t}) = 1 + \frac{q(\mathbf{t})}{r}, \qquad q(\mathbf{t}) =
(\mathbf{t}-\mu)^\top\Sigma^{-1}(\mathbf{t}-\mu), \tag{3} \\

\\ z(\mathbf{t}) = \sum\_{k=1}^{p}\frac{\alpha_k(t_k -
\mu_k)}{\sigma_k}, \tag{4} \\

\\ I(\mathbf{t}) = \int_0^\infty u^{(p+r)/2 - 1}\\e^{-u}\\
\Phi\\\left(z(\mathbf{t})\sqrt{\frac{2u}{A(\mathbf{t})\\r}}\right)du,
\tag{5} \\

and \\C(r,\Sigma)\\ collects the normalisation constants.

### 2.2 The Integral \\I(\mathbf{t})\\

The integral in Eq. (5) does not have a closed form for \\p \> 1\\ and
must be evaluated numerically. The integrand has the form
\\u^a\\e^{-u}\\\Phi(\cdot)\\, where \\a = (p+r)/2 - 1\\. The factor
\\u^a e^{-u}\\ is exactly the kernel of a Gamma distribution, which
makes **Gauss–Laguerre quadrature** the natural choice:

\\ I(\mathbf{t}) \approx \sum\_{q=1}^{Q}
w_q\\u_q^a\\\Phi\\\left(z(\mathbf{t})\sqrt{\frac{2u_q}{A(\mathbf{t})\\r}}\right),
\tag{6} \\

where \\(u_q, w_q)\_{q=1}^Q\\ are the standard Gauss–Laguerre nodes and
weights (i.e. the quadrature rule for \\\int_0^\infty
e^{-u}\\g(u)\\du\\), and \\Q = 32\\ in the implementation.

### 2.3 Why \\Q = 32\\?

The choice \\Q = 32\\ was validated by parity tests against R’s
`stats::integrate(rel.tol = 1e-12)` for \\p = 2\\ and \\\log r \in
\[\log 2, \log 100\]\\. For \\p = 2\\, the quadrature error is below
\\10^{-10}\\ across the tested range. For higher dimensions (\\p \geq
3\\), the exponent \\a\\ increases and additional nodes may be needed;
users should verify with the parity test suite.

### 2.4 Negative Log-Likelihood

\\ -\ell(\theta) = \frac{n}{2}\log\|\Sigma\| + \frac{np}{2}\log r +
n\\\log\Gamma\\\left(\frac{r}{2}\right) +
\frac{p+r}{2}\sum\_{i=1}^{n}\log A_i - \sum\_{i=1}^{n}\log I_i +
\text{penalty}, \tag{7} \\

where:

- \\\log\|\Sigma\| = 2\sum_k\log L\_{\text{cov}}(k,k)\\,
- \\A_i = 1 + q_i/r\\,
- \\I_i\\ is evaluated via Eq. (6),
- the penalty includes ridge terms on \\\log\sigma\\, \\\alpha\\, and
  \\\log r\\.

Note: the constants \\-(p/2 + 1)\\n\log 2\\ and
\\\frac{p}{2}\\n\log(2\pi)\\ are dropped to match the convention in
Models 1–4.

## 3 Parameterisation

The unconstrained parameter vector \\\theta\\ has length \\3p +
p(p-1)/2 + 1\\:

\\ \theta = \bigl\[\underbrace{\mu}\_{p},\\
\underbrace{\log\sigma}\_{p},\\
\underbrace{v\_{\text{C-vine}}}\_{p(p-1)/2},\\
\underbrace{\alpha}\_{p},\\ \underbrace{\log r}\_{1}\bigr\]. \\

| Block | Length | Constraint | Meaning |
|----|----|----|----|
| \\\mu\\ | \\p\\ | Unconstrained | Niche centroid |
| \\\log\sigma\\ | \\p\\ | \\\sigma_k \> 0\\ | Marginal standard deviations |
| \\v\\ | \\p(p-1)/2\\ | Unconstrained | C-vine partial correlations |
| \\\alpha\\ | \\p\\ | Unconstrained | Skewness parameters |
| \\\log r\\ | 1 | \\r \> 0\\ via \\\exp\\ | Degrees of freedom |

The degrees-of-freedom parameter \\r\\ is stored as \\\log r\\ at
position \\\theta\[3p + n_v\]\\ to enforce positivity.

## 4 Implementation Call Chain

    optimize_niche(likelihood = "skew_t")
      └─ niche_obj.cpp::create_niche_obj_ptr()
           └─ loglik_niche_math_skew_t_eigen()         [C++, lines 760–823]
                ├─ build_L_cov()                        [theta → mu, sigma, L_cov]
                ├─ alpha = Map(theta + 2*p + n_v, p)    [line 776]
                ├─ log_r = theta[3*p + n_v]; r = exp()  [lines 777–778]
                ├─ Mahalanobis q via forward solve       [lines 785–786]
                ├─ z_vec = skew_z(diff, sigma, alpha)    [line 789]
                ├─ log_det = 2 * sum log L_cov(k,k)     [lines 791–793]
                ├─ For each i:
                │    A_i = 1 + q_i/r                     [line 801]
                │    sum_log_A += log(A_i)               [line 802]
                │    sum_log_I += log_skew_t_integral()  [line 803]
                ├─ prior_penalty_value(has_alpha=true)   [lines 812–813]
                └─ neg_log = 0.5*n*log_det + 0.5*n*p*log_r
                             + n*lgamma(r/2)
                             + half_kpr*sum_log_A
                             - sum_log_I + penalty       [lines 814–819]

### 4.1 The Quadrature Engine (lines 725–752)

The function `log_skew_t_integral(q, z, r, k)` evaluates \\\log I\\ for
a single point:

1.  Computes \\A = 1 + q/r\\ and the exponent \\a = (p+r)/2 - 1\\ (lines
    727–728).
2.  Loops over \\Q = 32\\ Gauss–Laguerre nodes \\(u_q, w_q)\\ (lines
    733–744):
    - Argument to \\\Phi\\: \\\text{arg} = z\sqrt{2u_q/(Ar)}\\ (line
      740).
    - Log-term: \\\log w_q + a\log u_q + \log\Phi(\text{arg})\\ (line
      742).
3.  Applies the max-shift logsumexp trick across all \\Q\\ terms (lines
    746–751).

The nodes and weights are hard-coded as `constexpr double` arrays
(`kGLNodes`, `kGLWeights`, lines 677–715), generated by
`numpy.polynomial.laguerre.laggauss(32)`.

### 4.2 R Entry Point

``` r

optimize_niche(env_occ, likelihood = "skew_t",
               prior_log_sigma_lambda = 1.0,
               prior_alpha_lambda = 1.0)
```

## 5 Analytical Properties

| Property                                  | Value                           |
|-------------------------------------------|---------------------------------|
| Parameters                                | \\3p + p(p-1)/2 + 1\\           |
| Special case \\r \to \infty\\             | Skew-normal (Model 4)           |
| Special case \\\alpha = 0\\               | (Non-skewed) multivariate \\t\\ |
| Special case \\\alpha = 0, r \to \infty\\ | Multivariate normal (Model 1)   |
| Quadrature nodes                          | \\Q = 32\\ (Gauss–Laguerre)     |
| Background data                           | Not used                        |

## 6 When To Use

- When the species may have **heavy-tailed** environmental tolerances
  (survives extreme conditions more than a Gaussian predicts).
- As a **model-selection candidate**: compare AIC/BIC with
  `"skew_normal"` and `"presence_only"`. If \\\hat{r} \> 50\\, the
  skew-t is not adding value over the skew-normal.
- For species at **range margins** where extreme environmental events
  are ecologically relevant.

## 7 Limitations

1.  **Computational cost**. Each log-likelihood evaluation requires \\Q
    = 32\\ \\\Phi\\ evaluations per data point, making this model
    approximately \\32 \times\\ slower than the skew-normal per
    evaluation.

2.  **Quadrature accuracy for \\p \geq 3\\**. The parity tests validate
    \\Q = 32\\ for \\p = 2\\. For higher dimensions, the exponent \\a =
    (p+r)/2 - 1\\ increases, potentially requiring more nodes. Users
    should run the parity test suite (`test-loglik_math_cpp_parity.R`)
    for their specific \\p\\.

3.  **No background correction**. This is a presence-only model. For
    background-corrected skew-t, use `"skew_t_weighted"` (Model 7).

4.  **Parameter identifiability**. With \\3p + p(p-1)/2 + 1\\
    parameters, the model requires large \\n\\. For \\p = 2\\, this is 8
    parameters; for \\p = 5\\, it is 26 parameters.

## References

Branco, M. D., and D. K. Dey. 2001. “A General Class of Multivariate
Skew-Elliptical Distributions.” *Journal of Multivariate Analysis* 79
(1): 99–113. <https://doi.org/10.1006/jmva.2000.1960>.
