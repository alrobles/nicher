#' Validate a nicher object
#'
#' Internal function that checks whether an object has the expected structure
#' for a \code{nicher} object.  Raises an informative error if any required
#' component is missing or has the wrong type.
#'
#' @param x Object to validate.
#' @return \code{x} invisibly, if validation passes.
#' @keywords internal
check_nicher <- function(x) {
  if (!inherits(x, "nicher")) {
    stop("Object must be of class 'nicher'")
  }

  required <- c("loglik", "math_params", "bioscale_params",
                "optimx_table", "model", "method")
  missing_fields <- setdiff(required, names(x))
  if (length(missing_fields) > 0) {
    stop("nicher object is missing required fields: ",
         paste(missing_fields, collapse = ", "))
  }

  if (!is.numeric(x$loglik) || length(x$loglik) != 1) {
    stop("'loglik' must be a single numeric value")
  }

  if (!is.numeric(x$math_params)) {
    stop("'math_params' must be a numeric vector")
  }

  if (!is.list(x$bioscale_params)) {
    stop("'bioscale_params' must be a list")
  }

  bio_required <- c("mu", "L", "S", "variances")
  missing_bio <- setdiff(bio_required, names(x$bioscale_params))
  if (length(missing_bio) > 0) {
    stop("'bioscale_params' is missing fields: ",
         paste(missing_bio, collapse = ", "))
  }

  if (!is.character(x$model) || length(x$model) != 1) {
    stop("'model' must be a single character string")
  }

  valid_models <- c("unweighted", "weighted", "presence_only")
  if (!x$model %in% valid_models) {
    stop("'model' must be one of: ",
         paste(valid_models, collapse = ", "))
  }

  if (!is.character(x$method) || length(x$method) != 1) {
    stop("'method' must be a single character string")
  }

  invisible(x)
}
