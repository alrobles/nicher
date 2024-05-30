#' Ellipse constructor
#'
#' Constructor function that create an ellipse.
#' @param param_list A list of parameters of an ellipse object. First element
#' is the mu vector. Second element is the square matrix of covariates
#'
#' @return \code{ellipse}
#'
#' @export
#'
#' @examples
#' ellipse(get_ellip_par(spOccPnts))
#'
ellipse <- function(param_list) {

  # create list with the structure of a 'ellipse'.
  ellipse <- list(mu = param_list[[1]], S = param_list[[2]])

  # set class of object to "ellipse".
  class(ellipse) <- "ellipse"

  # return ellipse
  ellipse
}
