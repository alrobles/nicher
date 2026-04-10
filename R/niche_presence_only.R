#' Fit presence-only Gaussian niche model
#'
#' Fits a multivariate normal niche model using only presence points
#' (no background correction), using the compiled C++ backend (XPtr).
#'
#' Supports both:
#' \itemize{
#'   \item \strong{Single-start}: \code{start} is a numeric vector.
#'   \item \strong{Multi-start}: \code{start} is a list of numeric vectors;
#'         the best (lowest negative log-likelihood) is selected.
#' }
#'
#' @param occ Numeric matrix of environmental values at presence points.
#' @param eta LKJ prior parameter for the correlation structure (default 1).
#' @param start A numeric vector (single-start) or list of numeric vectors (multi-start).
#' @param ... Ignored. Present only for wrapper compatibility.
#'
#' @return A list with:
#'   \itemize{
#'     \item theta — best parameter vector (renamed from \code{par})
#'     \item value — negative log-likelihood
#'     \item conv — convergence code
#'     \item all_results — only present for multi-start
#'   }
#'
#' @examples
#' occ    <- as.matrix(example_env_occ_3d)
#' theta0 <- start_theta(example_env_occ_3d)
#' res    <- niche_presence_only(occ = occ, start = theta0)
#' res$value
#'
#' @export
niche_presence_only <- function(occ, eta = 1, start = NULL, ...) {
  # Build compiled C++ presence-only objective
  xptr <- create_niche_obj_ptr(
    env_occ    = as.matrix(occ),
    env_m      = NULL,
    likelihood = "presence_only",
    eta        = eta
  )

  # --- SINGLE START ---
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

  # --- MULTI-START ---
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

  stop("Start must be a numeric vector (single-start) or a list of numeric vectors (multi-start).")
}
