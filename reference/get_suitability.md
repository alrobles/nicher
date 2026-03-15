# get_suitability Function that calculates the log(suitability)

get_suitability Function that calculates the log(suitability)

## Usage

``` r
get_suitability(df, el_pars)
```

## Arguments

- df:

  A data frame with environmental information

- el_pars:

  A list with ellipse parameter. Contains a vector of means and a
  covariance matrix. Could came from a MLE after optimization or
  calculated from only presence data.

## Value

A vector with suitability

## Examples

``` r
pars <- get_negloglike_optimr_par(head(spOccPnts, 10), samMPts, fastmode = TRUE, itnmax = 1)
#> fastmode
#> TRUE
#> Optimization starting 
#> Methods used:
#> L-BFGS-B
#> Optimization finished. 
#> L-BFGS-B is the returned method
get_suitability(spOccPnts, pars)
#>  [1] 0.915758125 0.154292726 0.427894206 0.290579611 0.050198165 0.159878360
#>  [7] 0.061999401 0.416716329 0.291801927 0.577055231 0.025700023 0.058426207
#> [13] 0.244305004 0.407373661 0.048101978 0.923841938 0.034415976 0.029771727
#> [19] 0.026714557 0.071101181 0.046989809 0.009001979 0.039177911 0.013646957
#> [25] 0.010815686 0.031154608 0.019897935 0.006994209 0.007060055 0.045644363
#> [31] 0.024224954 0.020559529 0.565013414 0.116001391 0.090225498 0.290244925
#> [37] 0.231235693 0.216829972 0.021779331 0.014860112 0.021041879 0.120453224
#> [43] 0.107388617 0.014500360 0.019106568 0.714222740 0.025911605 0.069760227
#> [49] 0.040940437 0.220137542 0.266620629 0.053955070 0.117319370 0.136159929
#> [55] 0.377512734 0.038190510 0.267222654 0.290244925 0.838787608 0.190839532
#> [61] 0.197259791 0.565013414 0.116774376 0.037260418 0.923841938 0.023197909
#> [67] 0.047238927 0.072268249 0.033604180 0.052329153 0.092089718 0.703420330
#> [73] 0.841041997
```
