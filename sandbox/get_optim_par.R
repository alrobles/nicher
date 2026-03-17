#' get_optim_par Get the parameters for the optimization of
#' the Mahalanobis ellipse
#'
#' @param df A data frame with environmental information
#'
#' @return A list with two objects. A vector with centers of ellipse and a
#' matrix with the inverse of covariance matrix
#' @export
#'
#' @examples
#' get_optim_par(spOccPnts)
get_optim_par <- function(df){
  mu <- colMeans(df, na.rm = TRUE)
  A <-   get_A_matrix(df)
  par_list <- list(mu = mu, A = A)
  return(par_list)
}
