#' Print method for \code{nicher} objects
#'
#' Displays the key results of a fitted \code{nicher} model in a concise,
#' human-readable format using biological-scale parameters.  The full
#' mathematical-scale parameter vector and the optimizer comparison table are
#' accessible via \code{\link{summary.nicher}} or by direct list access (e.g.
#' \code{fit$math_params}).
#'
#' @param x   An object of class \code{"nicher"}.
#' @param ... Currently unused; present for S3 method consistency.
#'
#' @return Returns \code{x} invisibly.
#'
#' @seealso \code{\link{nicher}}, \code{\link{summary.nicher}}
#'
#' @export
#'
#' @examples
#' fit <- nicher(spOccPnts, samMPts, model = "presence_only")
#' print(fit)
print.nicher <- function(x, ...) {
  cat("Nicher model\n")
  cat("  Method :", x$method, "\n")
  cat("  Model  :", x$model,  "\n")
  cat("  Log-likelihood:", round(x$loglik, 4), "\n\n")

  cat("Optimum (mu):\n")
  print(round(x$bio_params$mu, 4))

  cat("\nStd. deviations (sigma):\n")
  print(round(x$bio_params$sigma, 4))

  cat("\nCorrelation matrix (R):\n")
  print(round(x$bio_params$R, 4))

  invisible(x)
}
