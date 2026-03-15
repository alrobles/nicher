#' Log-verosimilitud para modelo de nicho gaussiano (escala biológica)
#'
#' @param sam1 Data.frame o matriz de puntos de presencia (filas = observaciones, columnas = variables).
#' @param sam2 Data.frame o matriz de puntos de fondo ambiental.
#' @param mu Vector de medias (óptimo) de longitud p.
#' @param L Matriz triangular inferior de la descomposición de Cholesky de la matriz de covarianza S (S = L %*% t(L)).
#'          Las diagonales deben ser positivas.
#' @return Log-verosimilitud (escalar).
loglik_bio_niche <- function(sam1, sam2, mu, L) {
  # Convertir a matriz si es necesario
  sam1 <- as.matrix(sam1)
  sam2 <- as.matrix(sam2)
  p <- length(mu)

  # Calcular matriz de covarianza
  S <- L %*% t(L)
  logdet <- as.numeric(determinant(S, logarithm = TRUE)$modulus)

  # Función auxiliar para calcular distancias de Mahalanobis
  mahal <- function(x, mu, S) {
    # x es vector o matriz de una fila
    if (is.vector(x)) {
      x <- matrix(x, nrow = 1)
    }
    # Calcular (x - mu) %*% solve(S) %*% t(x - mu) para cada fila
    # Usamos la descomposición de Cholesky para eficiencia
    # Pero para simplicidad, usamos solve
    x_centered <- sweep(x, 2, mu)
    # Usamos la fórmula: diag(x_centered %*% solve(S) %*% t(x_centered))
    # Más eficiente: usar forwardsolve/backsolve con L
    # Pero aquí usamos directamente mahalanobis de stats
    stats::mahalanobis(x, center = mu, cov = S, inverted = FALSE)
  }

  q1 <- mahal(sam1, mu, S)
  q2 <- mahal(sam2, mu, S)

  n1 <- length(q1)
  n2 <- length(q2)

  # Log-verosimilitud: diferencia de log-densidades (sin constante 2pi)
  logL <- -0.5 * (n1 * logdet + sum(q1)) + 0.5 * (n2 * logdet + sum(q2))

  return(logL)
}
