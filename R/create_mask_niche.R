#' Create a parameter mask for the Gaussian niche model
#'
#' @param p Number of environmental variables.
#' @return Named numeric vector of NAs for mu and lower-triangular L parameters.
create_mask_niche <- function(p) {
  k <- p * (p + 1) / 2
  names_mu <- paste0("mu", 1:p)
  names_L <- paste0("L", 1:k)
  stats::setNames(rep(NA_real_, p + k), c(names_mu, names_L))
}
