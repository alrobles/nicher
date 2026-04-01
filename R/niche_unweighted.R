#' Fit unweighted Gaussian niche model
#'
#' Fits the classical unweighted ellipsoid niche model (Jiménez & Soberón 2019)
#' using the compiled C++ backend (XPtr). Supports single-start and multi-start
#' optimization.
#'
#' @param occ Numeric matrix of environmental values at presence points.
#' @param M Numeric matrix of background environmental values (same columns as occ).
#' @param eta LKJ prior parameter (default 1).
#' @param start Either a numeric vector (single start) or list of numeric vectors (multi-start).
#' @param ... Ignored.
#'
#' @return A list with:
#'   \itemize{
#'     \item theta — best parameter vector (renamed from \code{par})
#'     \item value — negative log-likelihood
#'     \item conv — convergence code
#'     \item all_results — present only for multi-start
#'   }
#'
#' @export
niche_unweighted <- function(occ, M, eta = 1, start = NULL, ...) {
  xptr <- create_niche_obj_ptr(
    env_occ    = as.matrix(occ),
    env_m      = as.matrix(M),
    likelihood = "unweighted",
    eta        = eta
  )

  # --- Single start ---
  if (is.numeric(start)) {
    res <- optimize_niche_xptr(
      start       = start,
      xptr        = xptr,
      multi_start = FALSE
    )
    res$theta <- res$par
    res$par <- NULL
    return(res)
  }

  # --- Multi-start ---
  if (is.list(start)) {
    res <- optimize_niche_xptr(
      start       = start,
      xptr        = xptr,
      multi_start = TRUE
    )
    res$theta <- res$par
    res$par <- NULL
    return(res)
  }

  stop("Start must be numeric (single-start) or list (multi-start).")
}
