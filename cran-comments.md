## New submission

This is the first CRAN release of `nicher`.

## Dependency note

`ucminfcpp` is now available on CRAN (>= 1.0.0), so `nicher` no longer
has a hard dependency-blocker on a non-CRAN package.

## Test environments

- local: Ubuntu 22.04 / 26.04, R 4.1.2 / 4.5.x
- GitHub Actions: ubuntu-latest (R-release, R-devel), macOS-latest
  (R-release), windows-latest (R-release)

## R CMD check results

1 error | 0 warnings | 2 notes

### ERROR

- `checking CRAN incoming feasibility ... ERROR`
  - `Conflicting package names (submitted: nicher, existing: nicheR)`
  - CRAN already has the package `nicheR`, whose name differs only by
    case from `nicher`. CRAN policy requires package names to be
    distinguishable, so `nicher` must be renamed before submission.

### NOTEs

- `checking installed package size ... NOTE`
  - installed size is 22.3Mb; `libs/` is 20.7Mb.
  - This is due to the compiled C++ code (Rcpp / RcppEigen / RcppParallel).

- `unable to verify current time` (local check only, not on CRAN).

### Possibly invalid URLs

`--as-cran` reports several DOI URLs returning 403. These are valid
DOIs and the 403 is a server response to the automated checker, not an
actual invalid URL.

## Outstanding blockers before CRAN submission

1. Rename the package so it no longer conflicts with the existing CRAN
   package `nicheR`.
2. Decide whether the installed size NOTE is acceptable or needs to be
   reduced (e.g. by stripping debug symbols or splitting functionality).
