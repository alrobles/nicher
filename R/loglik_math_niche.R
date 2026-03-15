#' Log-verosimilitud para modelo de nicho gaussiano en escala matemática
#'
#' @param param_vector Vector nombrado con parámetros en escala matemática.
#'                     Debe contener mu1..mup y L1..Lk en el orden de la triangular inferior por columnas.
#'                     Las diagonales de L deben estar en escala logarítmica.
#' @param sam1 Data.frame de puntos de presencia.
#' @param sam2 Data.frame de puntos de fondo.
#' @param mask Vector nombrado con valores fijos (opcional). Los nombres deben coincidir con los de param_vector.
#' @param negative Lógico; si TRUE devuelve el negativo de la log-verosimilitud.
#' @return Log-verosimilitud (o su negativo).
loglik_math_niche <- function(param_vector, sam1, sam2, mask = NULL, negative = TRUE) {
  # Determinar p a partir de los nombres de mu
  mu_names <- grep("^mu", names(param_vector), value = TRUE)
  if (length(mu_names) == 0) stop("No se encontraron parámetros mu en param_vector")
  p <- length(mu_names)

  # Número de elementos de L esperados
  k <- p * (p + 1) / 2
  nombres_completos <- c(paste0("mu", 1:p), paste0("L", 1:k))

  # Si hay máscara, completar el vector
  if (!is.null(mask)) {
    # Crear vector completo con NAs
    full <- stats::setNames(rep(NA_real_, length(nombres_completos)), nombres_completos)
    # Insertar valores de mask
    if (any(!names(mask) %in% nombres_completos)) {
      stop("Nombres en mask no válidos: ",
           paste(setdiff(names(mask), nombres_completos), collapse = ", "))
    }
    full[names(mask)] <- mask
    # Insertar valores de param_vector (que deben ser los libres)
    if (any(!names(param_vector) %in% nombres_completos)) {
      stop("Nombres en param_vector no válidos: ",
           paste(setdiff(names(param_vector), nombres_completos), collapse = ", "))
    }
    overlap <- intersect(names(param_vector), names(mask))
    if (length(overlap) > 0) {
      stop("Los nombres en param_vector y mask no pueden solaparse: ",
           paste(overlap, collapse = ", "))
    }
    full[names(param_vector)] <- param_vector
    if (anyNA(full)) {
      stop("Faltan parámetros: ", paste(names(full)[is.na(full)], collapse = ", "))
    }
    param_vector <- full
  }

  # Verificar que todos los nombres están presentes
  if (!all(nombres_completos %in% names(param_vector))) {
    stop("El vector de parámetros no contiene todos los nombres requeridos")
  }

  # Convertir a escala biológica
  bio <- math_to_bio_niche(param_vector)

  # Calcular log-verosimilitud
  logL <- loglik_bio_niche(sam1, sam2, bio$mu, bio$L)

  if (negative) return(-logL) else return(logL)
}
