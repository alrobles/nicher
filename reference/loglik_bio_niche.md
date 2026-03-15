# Log-verosimilitud para modelo de nicho gaussiano (escala biológica)

Log-verosimilitud para modelo de nicho gaussiano (escala biológica)

## Usage

``` r
loglik_bio_niche(sam1, sam2, mu, L)
```

## Arguments

- sam1:

  Data.frame o matriz de puntos de presencia (filas = observaciones,
  columnas = variables).

- sam2:

  Data.frame o matriz de puntos de fondo ambiental.

- mu:

  Vector de medias (óptimo) de longitud p.

- L:

  Matriz triangular inferior de la descomposición de Cholesky de la
  matriz de covarianza S (S = L Las diagonales deben ser positivas.

## Value

Log-verosimilitud (escalar).
