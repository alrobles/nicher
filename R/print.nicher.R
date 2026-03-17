#' Print method for nicher objects
#'
#' Prints a concise, human-readable summary of a fitted \code{nicher} object.
#' Shows the model type, optimization method, log-likelihood, and the
#' biological-scale optimized parameters (environmental optima and marginal
#' standard deviations).
#'
#' @param x An object of class \code{"nicher"}.
#' @param digits Integer; number of significant digits to print (default 4).
#' @param ... Further arguments passed to or from other methods (currently
#'   unused).
#' @return \code{x} invisibly.
#' @export
#' @method print nicher
#'
#' @examples
#' fit <- nicher(spOccPnts, samMPts, model = "presence_only")
#' print(fit)
print.nicher <- function(x, digits = 4, ...) {
  check_nicher(x)
  cat("Ecological niche model (nicher)\n")
  cat("  Model  :", x$model, "\n")
  cat("  Method :", x$method, "\n")
  cat("  Log-likelihood:", round(x$loglik, digits), "\n\n")

  cat("Biological-scale parameters:\n")

  mu        <- x$bioscale_params$mu
  variances <- x$bioscale_params$variances
  sds       <- sqrt(variances)

  env_names <- if (!is.null(x$env_names)) x$env_names else names(mu)
  if (is.null(env_names)) {
    env_names <- paste0("var", seq_along(mu))
  }

  tbl <- data.frame(
    variable = env_names,
    optimum  = round(mu, digits),
    sd       = round(sds, digits),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  print(tbl, row.names = FALSE)

  invisible(x)
}
