#' Negative log-likelihood (weighted-normal, math scale) with explicit index control
#'
#' This version extends the original \code{loglik_niche_math_weighted()} by allowing
#' the user to pass explicit subsampling indices (\code{den_idx}, \code{kde_idx})
#' and precomputed KDE weights (\code{precomp_w_den}) directly, enabling
#' high-performance workflows where costly KDE denominators are computed only once.
#'
#' Internally this function calls
#' \code{\link{loglik_niche_math_weighted_integrated}}, which routes the work to
#' the native C++ implementation \code{loglik_niche_weighted_integrated_cpp}.
#'
#' ## Parameter logic
#'
#' * If the user supplies \code{den_idx} or \code{kde_idx}, they are used as-is
#'   without generating any new subsampling indices.
#'
#' * If both \code{den_idx} and \code{kde_idx} are \code{NULL}, the function falls
#'   back to the original behavior: subsampling indices are created using
#'   \code{m_subsample} and \code{m_kde_subsample}, respectively.
#'
#' * If \code{precomp_w_den} is provided, the function avoids recomputing KDE for
#'   the denominator, enabling large speedups for the weighted model.
#'
#' @param theta Numeric vector of unconstrained parameters (math scale).
#' @param env_occ Data frame or matrix with environmental values at presence points.
#' @param env_m Data frame or matrix with environmental values from M.
#' @param eta Numeric, shape parameter of the LKJ prior for the correlation matrix.
#' @param neg Logical. If \code{TRUE} (default), returns negative log-likelihood.
#'
#' @param m_subsample Optional integer or fraction. If \code{den_idx} is not given,
#'   this defines the number (or fraction) of rows of \code{env_m} used for the
#'   denominator subsample.
#'
#' @param m_kde_subsample Optional integer or fraction. If \code{kde_idx} is not given,
#'   this defines the KDE reference subsample size (or fraction).
#'
#' @param seed Optional integer seed for reproducible subsampling.
#'
#' @param den_idx Optional integer vector of 1-based row indices for the denominator
#'   subset. If provided, no new denominator indices are generated.
#'
#' @param kde_idx Optional integer vector of 1-based row indices for the KDE reference
#'   subset. If provided, no new KDE indices are generated.
#'
#' @param precomp_w_den Optional numeric vector of precomputed denominator KDE weights.
#'   Must match the size of \code{den_idx}. If provided, KDE for the denominator is
#'   not recomputed.
#'
#' @param ... Additional arguments (ignored; provided for compatibility).
#'
#' @return A scalar numeric value containing the (negative) log-likelihood.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Using explicit indices and precomputed KDE weights for maximum speed
#' den_idx <- sample.int(nrow(example_env_m_2d), 2000)
#' kde_idx <- sample.int(nrow(example_env_m_2d), 5000)
#' pre_w <- kde_gaussian(
#'   example_env_m_2d[den_idx, ],
#'   example_env_m_2d[kde_idx, ]
#' )
#'
#' loglik_niche_math_weighted(
#'   theta = start_theta(example_env_occ_2d),
#'   env_occ = example_env_occ_2d,
#'   env_m = example_env_m_2d,
#'   den_idx = den_idx,
#'   kde_idx = kde_idx,
#'   precomp_w_den = pre_w
#' )
#' }
loglik_niche_math_weighted <- function(
  theta,
  env_occ,
  env_m,
  eta = 1,
  neg = TRUE,
  m_subsample = NULL,
  m_kde_subsample = NULL,
  seed = NULL,
  den_idx = NULL,
  kde_idx = NULL,
  precomp_w_den = NULL,
  ...
) {
  # Dimension checks
  p <- ncol(env_occ)
  if (p != ncol(env_m)) {
    stop("env_occ and env_m must have the same number of columns")
  }

  n_m <- nrow(env_m)

  # -------------------------------------------------------------------------
  # CASE 1: User provides explicit indices → use them directly
  # -------------------------------------------------------------------------
  if (!is.null(den_idx) || !is.null(kde_idx)) {
    return(
      loglik_niche_math_weighted_integrated(
        theta = theta,
        env_occ = env_occ,
        env_m = env_m,
        eta = eta,
        neg = neg,
        den_idx = den_idx,
        kde_idx = kde_idx,
        precomp_w_den = precomp_w_den
      )
    )
  }

  # -------------------------------------------------------------------------
  # CASE 2: No explicit indices → fall back to original subsampling behavior
  # -------------------------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)

  pick_size <- function(x, nmax) {
    if (is.null(x)) {
      return(NULL)
    }
    if (length(x) != 1L ||
      !is.numeric(x) ||
      !is.finite(x) ||
      x <= 0) {
      stop("m_subsample/m_kde_subsample must be a single positive number.")
    }

    if (x < 1) max(1L, floor(x * nmax)) else min(nmax, as.integer(round(x)))
  }

  n_den <- pick_size(m_subsample, n_m)
  den_idx <- if (!is.null(n_den) && n_den < n_m) sample.int(n_m, n_den) else NULL

  n_kde <- pick_size(m_kde_subsample, n_m)
  kde_idx <- if (!is.null(n_kde) && n_kde < n_m) sample.int(n_m, n_kde) else NULL

  # -------------------------------------------------------------------------
  # Delegate to integrated C++ version
  # -------------------------------------------------------------------------
  loglik_niche_math_weighted_integrated(
    theta         = theta,
    env_occ       = env_occ,
    env_m         = env_m,
    eta           = eta,
    neg           = neg,
    den_idx       = den_idx,
    kde_idx       = kde_idx,
    precomp_w_den = precomp_w_den
  )
}
