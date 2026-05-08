# Build Cholesky factor of a correlation matrix from C‑vine partial correlations

Constructs the lower‑triangular Cholesky factor \\\mathbf{L}\\ of a \\d
\times d\\ correlation matrix using the C‑vine method described in
Lewandowski, Kurowicka & Joe (2009). The factor satisfies \\\mathbf{R} =
\mathbf{L}\mathbf{L}^\top\\ with \\\operatorname{diag}(\mathbf{R}) =
1\\.

## Usage

``` r
cvine_cholesky(v, d, eta = 1)
```

## Arguments

- v:

  Numeric vector of unconstrained reals, one for each C‑vine edge. The
  length must be \\d(d-1)/2\\. For \\d = 2\\, a single value is
  required. The order follows a \*\*level‑major\*\* sequence: first all
  edges of level 1 (i.e., between variable 1 and each later variable),
  then edges of level 2 (between variable 2 and later variables,
  conditioned on variable 1), and so on.

- d:

  Integer, dimension of the target correlation matrix (\\d \ge 1\\).

- eta:

  Positive numeric shape parameter for the LKJ‑C‑vine prior (default
  `1`). \\\eta = 1\\ gives a uniform distribution over correlation
  matrices; larger values concentrate mass near the identity.

## Value

A \\d \times d\\ lower‑triangular matrix \\\mathbf{L}\\ with positive
diagonal entries such that \\\mathbf{L}\mathbf{L}^\top\\ is a valid
correlation matrix (unit diagonal, positive definite). For \\d = 1\\,
returns a \\1 \times 1\\ matrix with entry `1`.

## Details

The algorithm proceeds in three steps:

1.  Each element of `v` is mapped to the unit interval via the logistic
    (sigmoid) function, then transformed to a partial correlation on
    \\(-1, 1)\\ using the quantile function of a symmetric Beta
    distribution with shape \\\phi_k = \eta + (d - k - 2)/2\\ (for level
    \\k\\, 0‑indexed).

2.  The table of partial correlations is converted to unconditional
    correlations using the Yule–Kendall recursion (vine recursion).

3.  For each new row \\j\\ (starting from \\j = 2\\), the algorithm
    solves a triangular system to obtain the first \\j-1\\ entries of
    the row, and sets the diagonal entry to maintain unit row norm.

The resulting \\\mathbf{L}\\ can be used directly to construct the
correlation matrix (`tcrossprod(L)`), or to build a covariance matrix by
scaling with standard deviations.

## References

Lewandowski, D., Kurowicka, D., & Joe, H. (2009). Generating random
correlation matrices based on vines and extended onion method. \*Journal
of Multivariate Analysis\*, 100(9), 1989–2001.
[doi:10.1016/j.jmva.2009.04.008](https://doi.org/10.1016/j.jmva.2009.04.008)

## Examples

``` r
# For a 2x2 correlation matrix, we need one parameter
v <- 0.5
L <- cvine_cholesky(v, d = 2, eta = 1)
R <- tcrossprod(L)
print(R)
#>           [,1]      [,2]
#> [1,] 1.0000000 0.2449187
#> [2,] 0.2449187 1.0000000

# For 3 dimensions, we need 3 parameters (d*(d-1)/2 = 3)
v <- c(0.1, -0.2, 0.8)
L <- cvine_cholesky(v, d = 3)
R <- tcrossprod(L)
all.equal(diag(R), rep(1, 3)) # Should be TRUE
#> [1] TRUE
```
