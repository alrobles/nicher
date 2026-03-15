#' Crear una máscara para parámetros de modelo de nicho gaussiano
#'
#' @param p Número de variables ambientales.
#' @return Vector nombrado con NAs para los parámetros mu y los elementos de la triangular inferior de L.
create_mask_niche <- function(p) {
  k <- p * (p + 1) / 2
  names_mu <- paste0("mu", 1:p)
  names_L <- paste0("L", 1:k)
  stats::setNames(rep(NA_real_, p + k), c(names_mu, names_L))
}
