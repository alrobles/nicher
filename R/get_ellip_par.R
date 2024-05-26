#' Get ellipsoid parameters. A function to compute average and the inverse of covariance matrix from environmental data
#'
#' @param env A data frame containing environmental variables
#' @return A list with computed average of environmental variables and the covariance matrix
#' @export
#'
#' @examples
#' get_ellip_par(spOccPnts)
get_ellip_par <- function(env) {
  mu <- colMeans(env, na.rm = TRUE)
  S <- stats::cov(env)
  par_list <- list(mu = mu, S = S)
  return(par_list)
}

