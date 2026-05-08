# Reference (pure-R) negative log-likelihood for the NCST density

Evaluates the Non-Central Skew t (NCST) negative log-likelihood using
pure base-R arithmetic (no C++ calls except for the C-vine Cholesky).
This is the Stage 1 reference implementation following the SSDLC:
correctness is everything, speed is irrelevant.

## Usage

``` r
loglik_niche_math_ncst_r(theta, env_occ, eta = 1)
```

## Arguments

- theta:

  Numeric vector of length 3\*p + p\*(p-1)/2 + 1:
  `[xi(p), log_sigma(p), v(n_v), alpha(p), log_r(1)]`.

- env_occ:

  Numeric matrix (n x p) of occurrence data.

- eta:

  Scalar shape parameter for the C-vine prior (default 1).

## Value

Scalar negative log-likelihood (same constant-dropping convention as the
C++ kernel).

## Details

The NCST stochastic form (Hasan & Chen 2025, Definition 1) is: T = X /
sqrt(Y / r), X ~ SN_k(xi, Omega, alpha), Y ~ chi^2_r, where the location
xi enters INSIDE the skew-normal BEFORE chi-squared scaling. This
differs from Azzalini/Branco & Dey's skew-t where location is added
AFTER scaling.

## References

Hasan, M. R. & Chen, M.-H. (2025). Flexible Modeling of Multivariate
Skewed and Heavy-Tailed Data via a Non-Central Skew t Distribution.
arXiv:2507.10465v1.
