# Log-verosimilitud para modelo de nicho gaussiano en escala matemática

Log-verosimilitud para modelo de nicho gaussiano en escala matemática

## Usage

``` r
loglik_math_niche(param_vector, sam1, sam2, mask = NULL, negative = TRUE)
```

## Arguments

- param_vector:

  Vector nombrado con parámetros en escala matemática. Debe contener
  mu1..mup y L1..Lk en el orden de la triangular inferior por columnas.
  Las diagonales de L deben estar en escala logarítmica.

- sam1:

  Data.frame de puntos de presencia.

- sam2:

  Data.frame de puntos de fondo.

- mask:

  Vector nombrado con valores fijos (opcional). Los nombres deben
  coincidir con los de param_vector.

- negative:

  Lógico; si TRUE devuelve el negativo de la log-verosimilitud.

## Value

Log-verosimilitud (o su negativo).
