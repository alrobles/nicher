# R/optimize_niche_xptr.R
# Single-start niche optimization wrappers
#
# These functions provide a simple interface for single-start optimization
# of the three niche model likelihoods. Each function accepts a pre-specified
# starting parameter vector and returns the optimized log-likelihood together
# with the convergence code. They are designed for use in multi-start loops
# where starting points are generated externally (e.g., via
# start_theta_multiple() with the Sobol method).

# Default ucminf control parameters shared across all wrappers.
.default_ucminf_ctrl <- list(
  grad     = "central",
  gradstep = c(1e-6, 1e-8),
  grtol    = 1e-4,
  xtol     = 1e-8,
  stepmax  = 5,
  maxeval  = 2000
)

#' Single-start optimization: unweighted ellipsoid niche model
#'
#' Optimizes the unweighted log-likelihood (Jiménez & Soberón 2019) from a
#' single user-supplied starting vector using [ucminf::ucminf()].
#'
#' @param occ Numeric matrix (or data frame) of presence points
#'   (\code{n_occ x p}).
#' @param M Numeric matrix (or data frame) of background (accessibility area M)
#'   points (\code{n_m x p}).
#' @param start Numeric vector of starting parameters (math scale). Must have
#'   length \eqn{p + p + p(p-1)/2}.
#' @param eta Numeric. LKJ-C-vine prior shape parameter (default 1).
#' @param control List of control parameters passed to [ucminf::ucminf()].
#'   Missing entries are filled with package defaults.
#'
#' @return A list with:
#'   \item{value}{Maximized log-likelihood (positive, i.e., \code{-min}).}
#'   \item{conv}{Convergence code from [ucminf::ucminf()] (1 or 2 = success).}
#'   \item{theta}{Optimized parameter vector.}
#' @export
#'
#' @examples
#' \dontrun{
#' start <- start_theta(example_env_occ_3d)
#' res   <- niche_unweighted(
#'   occ   = as.matrix(example_env_occ_3d),
#'   M     = as.matrix(example_env_m_3d),
#'   start = start
#' )
#' cat("loglik =", res$value, "  conv =", res$conv, "\n")
#' }
niche_unweighted <- function(occ, M, start, eta = 1, control = list()) {
  occ   <- as.matrix(occ)
  M     <- as.matrix(M)
  start <- as.numeric(start)

  if (any(!is.finite(start))) {
    stop("All elements of 'start' must be finite.")
  }

  ctrl <- utils::modifyList(.default_ucminf_ctrl, control)

  fn <- function(theta) {
    loglik_niche_math_cpp(theta,
                          env_occ = occ,
                          env_m   = M,
                          eta     = eta,
                          neg     = TRUE)
  }

  out <- tryCatch(
    ucminf::ucminf(par     = start,
                   fn      = fn,
                   control = ctrl,
                   hessian = FALSE),
    error = function(e) {
      stop("Optimization failed in single-start mode: ", e$message)
    }
  )

  list(
    value = -as.numeric(out$value),
    conv  = as.integer(out$convergence),
    theta = out$par
  )
}

#' Single-start optimization: presence-only ellipsoid niche model
#'
#' Optimizes the presence-only log-likelihood (no background correction) from
#' a single user-supplied starting vector using [ucminf::ucminf()].
#'
#' @param occ Numeric matrix (or data frame) of presence points
#'   (\code{n_occ x p}).
#' @param start Numeric vector of starting parameters (math scale). Must have
#'   length \eqn{p + p + p(p-1)/2}.
#' @param eta Numeric. LKJ-C-vine prior shape parameter (default 1).
#' @param control List of control parameters passed to [ucminf::ucminf()].
#'   Missing entries are filled with package defaults.
#'
#' @return A list with:
#'   \item{value}{Maximized log-likelihood (positive, i.e., \code{-min}).}
#'   \item{conv}{Convergence code from [ucminf::ucminf()] (1 or 2 = success).}
#'   \item{theta}{Optimized parameter vector.}
#' @export
#'
#' @examples
#' \dontrun{
#' start <- start_theta(example_env_occ_3d)
#' res   <- niche_presence_only(
#'   occ   = as.matrix(example_env_occ_3d),
#'   start = start
#' )
#' cat("loglik =", res$value, "  conv =", res$conv, "\n")
#' }
niche_presence_only <- function(occ, start, eta = 1, control = list()) {
  occ   <- as.matrix(occ)
  start <- as.numeric(start)

  if (any(!is.finite(start))) {
    stop("All elements of 'start' must be finite.")
  }

  ctrl <- utils::modifyList(.default_ucminf_ctrl, control)

  fn <- function(theta) {
    loglik_niche_math_presence_only(theta,
                                    env_occ = occ,
                                    eta     = eta,
                                    neg     = TRUE)
  }

  out <- tryCatch(
    ucminf::ucminf(par     = start,
                   fn      = fn,
                   control = ctrl,
                   hessian = FALSE),
    error = function(e) {
      stop("Optimization failed in single-start mode: ", e$message)
    }
  )

  list(
    value = -as.numeric(out$value),
    conv  = as.integer(out$convergence),
    theta = out$par
  )
}

#' Single-start optimization: weighted ellipsoid niche model
#'
#' Optimizes the weighted log-likelihood (Jiménez & Soberón 2022) from a
#' single user-supplied starting vector using [ucminf::ucminf()]. Precomputed
#' KDE denominator weights can be supplied to avoid recomputing them on every
#' function evaluation, which significantly reduces runtime in multi-start
#' scenarios.
#'
#' @param occ Numeric matrix (or data frame) of presence points
#'   (\code{n_occ x p}).
#' @param M Numeric matrix (or data frame) of background (accessibility area M)
#'   points (\code{n_m x p}).
#' @param den_idx Integer vector of 1-based row indices into \code{M} selecting
#'   the denominator subsample.
#' @param kde_idx Integer vector of 1-based row indices into \code{M} selecting
#'   the KDE reference subsample.
#' @param precomp_w_den Numeric vector of precomputed KDE weights for the
#'   denominator points (i.e., \code{kde_gaussian(M[den_idx, ], M[kde_idx, ])}).
#'   Must have the same length as \code{den_idx}.
#' @param start Numeric vector of starting parameters (math scale). Must have
#'   length \eqn{p + p + p(p-1)/2}.
#' @param eta Numeric. LKJ-C-vine prior shape parameter (default 1).
#' @param control List of control parameters passed to [ucminf::ucminf()].
#'   Missing entries are filled with package defaults.
#'
#' @return A list with:
#'   \item{value}{Maximized log-likelihood (positive, i.e., \code{-min}).}
#'   \item{conv}{Convergence code from [ucminf::ucminf()] (1 or 2 = success).}
#'   \item{theta}{Optimized parameter vector.}
#' @export
#'
#' @examples
#' \dontrun{
#' occ3     <- as.matrix(example_env_occ_3d)
#' M3       <- as.matrix(example_env_m_3d)
#' set.seed(123)
#' den_idx  <- sample.int(nrow(M3), 3000)
#' kde_idx  <- sample.int(nrow(M3), 6000)
#' pre      <- kde_gaussian(M3[den_idx, ], M3[kde_idx, ])
#' start    <- start_theta(example_env_occ_3d)
#' res      <- niche_weighted(
#'   occ           = occ3,
#'   M             = M3,
#'   den_idx       = den_idx,
#'   kde_idx       = kde_idx,
#'   precomp_w_den = pre,
#'   start         = start
#' )
#' cat("loglik =", res$value, "  conv =", res$conv, "\n")
#' }
niche_weighted <- function(occ, M, den_idx, kde_idx, precomp_w_den,
                           start, eta = 1, control = list()) {
  occ           <- as.matrix(occ)
  M             <- as.matrix(M)
  den_idx       <- as.integer(den_idx)
  kde_idx       <- as.integer(kde_idx)
  precomp_w_den <- as.numeric(precomp_w_den)
  start         <- as.numeric(start)

  if (any(!is.finite(start))) {
    stop("All elements of 'start' must be finite.")
  }
  if (length(precomp_w_den) != length(den_idx)) {
    stop("'precomp_w_den' must have the same length as 'den_idx'.")
  }

  ctrl <- utils::modifyList(.default_ucminf_ctrl, control)

  fn <- function(theta) {
    loglik_niche_math_weighted_integrated(
      theta,
      env_occ       = occ,
      env_m         = M,
      eta           = eta,
      neg           = TRUE,
      den_idx       = den_idx,
      kde_idx       = kde_idx,
      precomp_w_den = precomp_w_den
    )
  }

  out <- tryCatch(
    ucminf::ucminf(par     = start,
                   fn      = fn,
                   control = ctrl,
                   hessian = FALSE),
    error = function(e) {
      stop("Optimization failed in single-start mode: ", e$message)
    }
  )

  list(
    value = -as.numeric(out$value),
    conv  = as.integer(out$convergence),
    theta = out$par
  )
}
