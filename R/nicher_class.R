# R/nicher_class.R
# Minimal S3 class system for nicher optimization results.
#
# Exports: new_nicher, print.nicher, assess, assess.nicher

# -----------------------------------------------------------------------
# Constructor
# -----------------------------------------------------------------------

#' Create a nicher optimization result object
#'
#' Constructs an S3 object of class \code{"nicher"} from the outputs of
#' \code{\link{optimize_niche}}.
#'
#' @param solutions A data frame with columns \code{start_id},
#'   \code{loglik}, \code{convergence}, and \code{full_par} (list column).
#' @param best A list with the best solution:
#'   \code{theta}, \code{loglik}, \code{convergence}.
#' @param likelihood Character. One of \code{"weighted"} or
#'   \code{"presence_only"}.
#' @param n_starts Integer. Total number of starting points used.
#' @param var_names Optional character vector of length \code{p} naming
#'   the environmental variables in the order assumed by
#'   \code{best$theta}. Stored on the object so
#'   \code{\link{predict.nicher}} can reorder a future
#'   \code{\link[terra]{SpatRaster}} by layer name.
#' @param eta Numeric scalar. The Beta-shape parameter passed to
#'   \code{\link{cvine_cholesky}} during fitting (default \code{1.0}).
#'   Stored on the object so \code{\link{predict.nicher}} reconstructs
#'   the same correlation matrix from \code{best$theta} that the
#'   optimizer used. Defaults to \code{1.0} (the LKJ-uniform prior),
#'   matching \code{\link{optimize_niche}}.
#'
#' @return An object of class \code{"nicher"}.
#' @keywords internal
new_nicher <- function(solutions, best, likelihood, n_starts,
                       var_names = NULL, eta = 1.0) {
  if (!is.numeric(eta) || length(eta) != 1L || !is.finite(eta) || eta <= 0) {
    stop("`eta` must be a single positive finite number.")
  }
  structure(
    list(
      solutions  = solutions,
      best       = best,
      likelihood = likelihood,
      n_starts   = as.integer(n_starts),
      var_names  = if (is.null(var_names)) NULL else as.character(var_names),
      eta        = as.numeric(eta)
    ),
    class = "nicher"
  )
}

# -----------------------------------------------------------------------
# print method
# -----------------------------------------------------------------------

#' Print a nicher object
#'
#' Displays a concise summary of a \code{"nicher"} optimization result,
#' including the likelihood type, number of starts, convergence count,
#' and the best log-likelihood achieved.
#'
#' @param x A \code{"nicher"} object returned by \code{\link{optimize_niche}}.
#' @param ... Ignored.
#'
#' @return Invisibly returns \code{x}.
#' @export
#' @examples
#' \dontrun{
#' result <- optimize_niche(
#'   env_occ    = example_env_occ_2d,
#'   env_m      = example_env_m_2d,
#'   num_starts = 10L,
#'   likelihood = "weighted"
#' )
#' print(result)
#' }
print.nicher <- function(x, ...) {
  n_conv <- sum(
    x$solutions$convergence %in% c(1L, 2L),
    na.rm = TRUE
  )
  cat("-- nicher optimization result --\n")
  cat("  Likelihood :", x$likelihood, "\n")
  cat(
    "  Starts     :", x$n_starts,
    "(converged:", n_conv, ")\n"
  )
  cat("  Best loglik:", round(x$best$loglik, 6L), "\n")
  cat("  Convergence:", x$best$convergence, "\n")
  if (!is.null(x$eta)) {
    cat("  eta        :", x$eta, "\n")
  }
  if (!is.null(x$var_names)) {
    cat("  Variables  :", paste(x$var_names, collapse = ", "), "\n")
  }
  invisible(x)
}

# -----------------------------------------------------------------------
# assess generic
# -----------------------------------------------------------------------

#' Assess acceptance criteria for an optimization result
#'
#' Generic function for evaluating whether an optimization result meets
#' acceptance criteria for niche model quality.
#'
#' @param x An object for which an \code{assess} method is defined.
#' @param ... Additional arguments passed to methods.
#'
#' @return A list of diagnostics and acceptance flags. See
#'   \code{\link{assess.nicher}} for details.
#' @export
assess <- function(x, ...) UseMethod("assess")

# -----------------------------------------------------------------------
# assess.nicher
# -----------------------------------------------------------------------

#' Assess acceptance criteria for a nicher result
#'
#' Evaluates the quality of a multi-start optimization by comparing
#' converged solutions. Returns one of four diagnostic flags:
#' \describe{
#'   \item{\code{"accepted_global"}}{Multiple starts agree closely on the
#'     same high log-likelihood value — the global optimum is likely found.}
#'   \item{\code{"accepted_noise"}}{Converged solutions are close but show
#'     minor numerical noise — likely a good solution.}
#'   \item{\code{"suggest_average"}}{Converged solutions spread across
#'     distinct values — consider averaging \code{theta} or increasing
#'     \code{num_starts}.}
#'   \item{\code{"needs_more_starts"}}{Too few converged solutions to draw
#'     conclusions — increase \code{num_starts}.}
#' }
#'
#' @param x A \code{"nicher"} object returned by \code{\link{optimize_niche}}.
#' @param tol_gap Numeric. Relative tolerance used to decide between
#'   \code{"accepted_global"} and \code{"accepted_noise"}. Default 0.01.
#' @param tol_dist Numeric. Relative tolerance used to decide between
#'   \code{"accepted_noise"} and \code{"suggest_average"}. Default 0.05.
#' @param min_converged Integer. Minimum number of converged solutions
#'   required to avoid \code{"needs_more_starts"}. Default 2.
#' @param ... Ignored.
#'
#' @return A named list with:
#' \describe{
#'   \item{\code{flag}}{Character scalar, one of the four flags above.}
#'   \item{\code{recommendation}}{Human-readable action recommendation.}
#'   \item{\code{gap}}{Absolute log-likelihood gap between the best and
#'     second-best converged solutions.}
#'   \item{\code{rel_gap}}{Relative gap (\code{gap / |best_loglik|}).}
#'   \item{\code{n_converged}}{Number of converged solutions.}
#'   \item{\code{best_loglik}}{Best log-likelihood value.}
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' result <- optimize_niche(
#'   env_occ    = example_env_occ_2d,
#'   env_m      = example_env_m_2d,
#'   num_starts = 20L,
#'   likelihood = "weighted"
#' )
#' diag <- assess(result)
#' cat("Flag          :", diag$flag, "\n")
#' cat("Recommendation:", diag$recommendation, "\n")
#' cat("Best log-lik  :", round(diag$best_loglik, 4L), "\n")
#' }
assess.nicher <- function(x,
                          tol_gap      = 0.01,
                          tol_dist     = 0.05,
                          min_converged = 2L,
                          ...) {
  if (min_converged < 2L) {
    stop("`min_converged` must be at least 2 because `gap` is computed from the top two converged solutions.")
  }
  sols   <- x$solutions

  # Filter to converged solutions (ucminf codes 1 or 2 = success)
  conv   <- sols[sols$convergence %in% c(1L, 2L), ]
  n_conv <- nrow(conv)

  # Derive best_loglik from the converged set; fall back to overall best
  # only when no solution converged (so the needs_more_starts path can
  # still report a meaningful value).
  best_ll <- if (n_conv > 0L) max(conv$loglik) else x$best$loglik

  if (n_conv < min_converged) {
    return(list(
      flag           = "needs_more_starts",
      recommendation = paste0(
        "Fewer than ", min_converged, " solutions converged. ",
        "Increase num_starts."
      ),
      gap         = NA_real_,
      rel_gap     = NA_real_,
      n_converged = n_conv,
      best_loglik = best_ll
    ))
  }

  # Gap: difference between best and second-best converged log-lik
  sorted_ll <- sort(conv$loglik, decreasing = TRUE)
  gap <- sorted_ll[1L] - sorted_ll[2L]

  # Relative gap (guard against best_ll = 0)
  rel_gap <- abs(gap) / (abs(best_ll) + .Machine$double.eps)

  if (rel_gap <= tol_gap) {
    flag <- "accepted_global"
    rec  <- paste0(
      "Global optimum likely found: ",
      "multiple converged starts agree closely."
    )
  } else if (rel_gap <= tol_dist) {
    flag <- "accepted_noise"
    rec  <- paste0(
      "Good solution found with minor numerical noise ",
      "across starts."
    )
  } else {
    flag <- "suggest_average"
    rec  <- paste0(
      "Converged solutions spread widely. ",
      "Consider averaging theta vectors or increasing num_starts."
    )
  }

  list(
    flag           = flag,
    recommendation = rec,
    gap            = gap,
    rel_gap        = rel_gap,
    n_converged    = n_conv,
    best_loglik    = best_ll
  )
}
