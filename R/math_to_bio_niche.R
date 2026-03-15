#' Convertir parámetros de escala matemática a biológica para modelo de nicho
#'
#' @param param_vector Vector nombrado con parámetros en escala matemática.
#'                     Nombres: mu1..mup, L1..Lk (k = p*(p+1)/2).
#'                     Las diagonales de L deben estar en escala logarítmica.
#' @return Lista con elementos `mu` (vector) y `L` (matriz triangular inferior).
math_to_bio_niche <- function(param_vector) {
  # Extraer mu
  mu <- param_vector[grep("^mu", names(param_vector))]
  p <- length(mu)

  # Extraer elementos de L
  L_elements <- param_vector[grep("^L", names(param_vector))]
  k <- length(L_elements)
  if (k != p * (p + 1) / 2) {
    stop("Número incorrecto de elementos para L")
  }

  # Reconstruir L en orden por columnas
  L <- matrix(0, p, p)
  L[lower.tri(L, diag = TRUE)] <- L_elements
  # Aplicar exponencial a las diagonales
  diag(L) <- exp(diag(L))

  list(mu = mu, L = L)
}
