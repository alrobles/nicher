#' Summary method for nicher objects
#'
#' Prints a detailed summary of a fitted \code{nicher} object, including
#' the model type, maximized log-likelihood, the full covariance matrix in
#' biological scale, and the complete optimization results table.
#'
#' @param object An object of class \code{"nicher"}.
#' @param digits Integer; number of significant digits (default 4).
#' @param ... Further arguments passed to or from other methods (currently
#'   unused).
#' @return An object of class \code{"summary.nicher"} (a list) containing
#'   the processed summary information, returned invisibly.
#' @export
#' @method summary nicher
#'
#' @examples
#' fit <- nicher(spOccPnts, samMPts, model = "presence_only")
#' summary(fit)
summary.nicher <- function(object, digits = 4, ...) {
  check_nicher(object)

  cat("=== nicher model summary ===\n\n")
  cat("Model   :", object$model, "\n")
  cat("Method  :", object$method, "\n")
  cat("Log-likelihood:", round(object$loglik, digits), "\n\n")

  # --- Environmental optima ---
  cat("Environmental optima (mu):\n")
  mu_named <- object$bioscale_params$mu
  env_names <- if (!is.null(object$env_names)) object$env_names else names(mu_named)
  if (is.null(env_names)) env_names <- paste0("var", seq_along(mu_named))
  names(mu_named) <- env_names
  print(round(mu_named, digits))
  cat("\n")

  # --- Covariance matrix ---
  cat("Covariance matrix (S):\n")
  S <- object$bioscale_params$S
  rownames(S) <- env_names
  colnames(S) <- env_names
  print(round(S, digits))
  cat("\n")

  # --- Correlation matrix ---
  cat("Correlation matrix (R):\n")
  R <- object$bioscale_params$R
  rownames(R) <- env_names
  colnames(R) <- env_names
  print(round(R, digits))
  cat("\n")

  # --- Math-scale parameters ---
  cat("Math-scale parameters (for reference):\n")
  print(round(object$math_params, digits))
  cat("\n")

  # --- Optimization table ---
  cat("Optimization table (all methods tried):\n")
  p <- length(object$bioscale_params$mu)
  k <- p * (p + 1) / 2
  n_par <- p + k
  extra_cols <- setdiff(names(object$optimx_table), names(object$math_params))
  tbl_display <- object$optimx_table[, extra_cols, drop = FALSE]
  print(tbl_display)

  out <- list(
    model          = object$model,
    method         = object$method,
    loglik         = object$loglik,
    mu             = object$bioscale_params$mu,
    S              = object$bioscale_params$S,
    R              = R,
    math_params    = object$math_params,
    optimx_table   = object$optimx_table
  )
  class(out) <- "summary.nicher"
  invisible(out)
}
