# Resolve KDE subsampling indices and precompute KDE weights ONCE per fit.

Implements the package's KDE sampling policy:

- Hard cap of 10000 distinct background combinations.

- Heuristic minimum representative sample of `max(500, 50 * 2^p)` rows;
  below this a [`warning()`](https://rdrr.io/r/base/warning.html) is
  fired.

- Optional fraction (`x < 1`) or integer count for `m_subsample` /
  `m_kde_subsample`.

## Usage

``` r
.resolve_weighted_inputs(
  env_occ,
  env_m,
  m_subsample = NULL,
  m_kde_subsample = NULL,
  seed = NULL
)
```

## Details

Returns a list with: `den_idx`, `kde_idx`, `w_occ`, `w_den`, `n_m`.
