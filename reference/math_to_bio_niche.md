# Convertir parámetros de escala matemática a biológica para modelo de nicho

Convertir parámetros de escala matemática a biológica para modelo de
nicho

## Usage

``` r
math_to_bio_niche(param_vector)
```

## Arguments

- param_vector:

  Vector nombrado con parámetros en escala matemática. Nombres:
  mu1..mup, L1..Lk (k = p\*(p+1)/2). Las diagonales de L deben estar en
  escala logarítmica.

## Value

Lista con elementos \`mu\` (vector) y \`L\` (matriz triangular
inferior).
