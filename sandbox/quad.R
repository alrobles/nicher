#' quad  Function that calculates quadratic terms
#'
#' @param xi An environmental data frame
#' @param mu A vector with the average of environmental variables
#' @param A matrix of inverse of covariance of environmental variables
#'
#' @return This function return a matrix with the quadratic form
#'
quad <- function(xi, mu, A) {
  if(is.null(dim(xi))){
    ax = xi - mu
  } else {
    ax <- apply(xi, 1, function(x) x - mu)
  }
  t(ax) %*% A %*% ax
}
