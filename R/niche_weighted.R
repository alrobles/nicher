# R/niche_weighted.R (updated)

#' Fit weighted-normal niche model (Jiménez & Soberón 2022)
#'
#' This function fits the *weighted-normal* niche model using a compiled
#' C++ backend via an external pointer (`XPtr`) created by
#' \code{\link{create_niche_obj_ptr}}.
#'
#' The weighted-normal model corrects for environmental bias in the
#' accessibility area \eqn{M} by incorporating kernel-density weights (KDE).
#'
#' For maximum performance and numerical stability, the KDE denominators
#' must be *precomputed* by the user and passed through
#' \code{precomp_w_den}. No KDE is recomputed inside the optimizer.
#'
#' This wrapper supports both:
#' \itemize{
#'   \item \strong{Single-start}: \code{start} is a numeric vector.
#'   \item \strong{Multi-start}: \code{start} is a list of numeric vectors.
#'         The best result (lowest negative log-likelihood) is returned.
#' }
#'
#' @param occ Numeric matrix of environmental values at presence points
#'   (rows = occurrences, columns = environmental variables).
#' @param M Numeric matrix of environmental values sampled from the
#'   accessibility area \eqn{M}. Must have the same number of columns as \code{occ}.
#' @param den_idx Integer vector of 1-based row indices selecting the
#'   denominator subset \eqn{M_\mathrm{den}}. Must have the same length as
#'   \code{precomp_w_den}.
#' @param kde_idx Integer vector of 1-based row indices selecting the KDE
#'   reference subset \eqn{M_\mathrm{kde}}. Must have length >= number of columns.
#' @param precomp_w_den Numeric vector of precomputed KDE weights matching
#'   \code{den_idx} in length.
#' @param eta Numeric LKJ shape parameter for the C-vine prior (default = 1).
#' @param start Starting value(s) for optimization:
#'   \itemize{
#'     \item numeric vector → single-start
#'     \item list of numeric vectors → multi-start
#'   }
#' @param ... Ignored. Present only for compatibility with other wrappers.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{theta} — best parameter vector (renamed from \code{par} for consistency with earlier versions)
#'   \item \code{value} — negative log-likelihood at optimum
#'   \item \code{conv} — convergence code
#'   \item \code{all_results} — data frame of all starts (multi-start only)
#' }
#'
#' @export
niche_weighted <- function(occ, M, den_idx, kde_idx, precomp_w_den,
                           eta = 1, start = NULL, ...) {
  # Validate precomputed denominators
  stopifnot(length(den_idx) == length(precomp_w_den))

  # Create compiled C++ objective function
  xptr <- create_niche_obj_ptr(
    env_occ        = as.matrix(occ),
    env_m          = as.matrix(M),
    likelihood     = "weighted",
    den_idx        = den_idx,
    kde_idx        = kde_idx,
    precomp_w_den  = precomp_w_den,
    eta            = eta
  )

  # --- SINGLE START (numeric vector) ---
  if (is.numeric(start)) {
    res <- optimize_niche_xptr(
      start       = start,
      xptr        = xptr,
      multi_start = FALSE
    )
    # Rename $par to $theta for consistency with previous versions
    res$theta <- res$par
    res$par <- NULL
    return(res)
  }

  # --- MULTI-START (list of numeric vectors) ---
  if (is.list(start)) {
    res <- optimize_niche_xptr(
      start       = start,
      xptr        = xptr,
      multi_start = TRUE
    )
    # For multi-start, $par is already the best; rename to $theta
    res$theta <- res$par
    res$par <- NULL
    return(res)
  }

  stop("Start must be a numeric vector (single-start) or a list (multi-start).")
}
