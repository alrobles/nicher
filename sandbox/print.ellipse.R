#' Print method for "ellipse" class
#'
#' @param x An object of class "ellipse"
#' @param ... Other arguments passed to or from other methods
#'
#' @export print.ellipse
#' @export
print.ellipse <- function(x, ...) {
  print(unclass(x))
}
