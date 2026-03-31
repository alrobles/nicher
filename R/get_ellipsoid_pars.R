#' Get ellipsoid parameters. A function to compute average and the inverse of covariance matrix from environmental data
#'
#' @param env A data frame containing environmental variables
#' @return A list with computed average of environmental variables and the
#' covariance matrix
#' @export
#'
#' @examples
#' get_ellipsoid_pars(example_env_occ_2d)
get_ellipsoid_pars <- function(env) {
  mu <- colMeans(env, na.rm = TRUE)
  s_mat <- stats::cov(env)
  list(mu = mu, s_mat = s_mat)
}
