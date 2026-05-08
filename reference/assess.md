# Assess acceptance criteria for an optimization result

Generic function for evaluating whether an optimization result meets
acceptance criteria for niche model quality.

## Usage

``` r
assess(x, ...)
```

## Arguments

- x:

  An object for which an `assess` method is defined.

- ...:

  Additional arguments passed to methods.

## Value

A list of diagnostics and acceptance flags. See
[`assess.nicher`](https://alrobles.github.io/nicher/reference/assess.nicher.md)
for details.
