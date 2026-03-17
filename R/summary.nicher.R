#' Summary method for \code{nicher} objects
#'
#' Prints a comprehensive summary of a fitted \code{nicher} model, including
#' biological-scale parameters (niche optimum, covariance matrix, standard
#' deviations, correlation matrix), the full mathematical-scale parameter
#' vector, and the optimizer comparison table.
#'
#' @param object An object of class \code{"nicher"}.
#' @param ...    Currently unused; present for S3 method consistency.
#'
#' @return Returns \code{object} invisibly.
#'
#' @seealso \code{\link{nicher}}, \code{\link{print.nicher}}
#'
#' @export
#'
#' @examples
#' fit <- nicher(spOccPnts, samMPts, model = "presence_only")
#' summary(fit)
summary.nicher <- function(object, ...) {
  cat("=== nicher model summary ===\n\n")
  cat("Method         :", object$method, "\n")
  cat("Model          :", object$model,  "\n")
  cat("Log-likelihood :", round(object$loglik, 6), "\n\n")

  cat("--- Biological-scale parameters ---\n\n")

  cat("Optimum (mu):\n")
  print(round(object$bio_params$mu, 6))

  cat("\nCovariance matrix (S):\n")
  print(round(object$bio_params$S, 6))

  cat("\nStd. deviations (sigma):\n")
  print(round(object$bio_params$sigma, 6))

  cat("\nCorrelation matrix (R):\n")
  print(round(object$bio_params$R, 6))

  cat("\n--- Mathematical-scale parameters ---\n\n")
  print(round(object$math_params, 6))

  cat("\n--- Optimization table ---\n\n")
  print(object$optim_table)

  invisible(object)
}
