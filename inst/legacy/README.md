# Legacy Functions — nicher

This folder archives R functions that have been removed from the active package
API as part of the `cleanup-agent` cleanup (roadmap Phase 1).  They are kept
here **for reference and backward-compatibility documentation only** and are
**not loaded** when the package is installed.

---

## `negloglike_multivariable.R`

**Status:** Deprecated — replaced by `loglik_presenceonly_math()`

**Purpose:**  
A backward-compatible wrapper for the presence-only negative log-likelihood.
It accepted `(mu, S, sam1, sam2)` and internally delegated to
`loglik_presenceonly_math(sam1, sam2, mu, S)`.

**Migration:**
```r
# Old (deprecated)
negloglike_multivariable(mu, S, sam1, sam2)

# New
loglik_presenceonly_math(sam1, sam2, mu, S)
```

---

## `get_optim_par.R`

**Status:** Dead code — superseded by `get_ellip_par()`

**Purpose:**  
Returned ellipsoid parameters in the format `list(mu, A)` where `A` is the
*inverse* covariance matrix (precision matrix).  Used exclusively with the
now-removed `get_negative_log()` function.

**Migration:**
```r
# Old (dead code)
par <- get_optim_par(df)   # returns list(mu, A) — A is inverse covariance

# New
par <- get_ellip_par(df)   # returns list(mu, S) — S is covariance matrix
```

---

## `get_negative_log.R`

**Status:** Dead code — superseded by `loglik_presenceonly_math()` and
`loglik_presenceonly_cpp()`

**Purpose:**  
Computed the negative log-likelihood directly from pre-computed Mahalanobis
quadratic forms `q1` (presence) and `q2` (background):

```r
get_negative_log(q1, q2)
# = 0.5 * sum(q1) + n * log(sum(exp(-0.5 * q2)))
```

Required callers to pre-compute `q1`/`q2` using `mahalanobis(..., inverted = TRUE)`
with the *inverse* covariance matrix from `get_optim_par()`.

**Migration:**
```r
# Old (dead code)
par <- get_optim_par(df)
q1  <- mahalanobis(sam1, par$mu, par$A, inverted = TRUE)
q2  <- mahalanobis(sam2, par$mu, par$A, inverted = TRUE)
get_negative_log(q1, q2)

# New
par <- get_ellip_par(df)
loglik_presenceonly_math(sam1, sam2, par$mu, par$S)
# or, for C++ performance:
loglik_presenceonly_cpp(as.matrix(sam1), as.matrix(sam2), par$mu, par$S)
```
