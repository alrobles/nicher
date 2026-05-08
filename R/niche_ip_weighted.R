#' Fit inverse-probability-weighted (IPW) normal niche model
#'
#' Fits the inverse-probability-weighted (IPW) normal niche model using a
#' compiled C++ backend via an external pointer for fast evaluation.
#'
#' This model assumes the data-generating process of Jimenez et al. (2019):
#' presence records are drawn from the fundamental niche density
#' \eqn{f(\mathbf{x};\mu,\Sigma)} restricted to the environments available
#' in \eqn{M}.  Because background grid cells are not uniformly distributed
#' in environmental space (their density is \eqn{g(\mathbf{x})}), a naive
#' sum over the background would estimate
#' \eqn{\int f\,g\,d\mathbf{x}} rather than the required
#' \eqn{\int f\,d\mathbf{x}}.  The IPW correction divides each background
#' term by \eqn{\hat{g}} (a kernel density estimate of the environmental
#' density of \eqn{M}), yielding a Horvitz--Thompson estimator
#' (Horvitz & Thompson 1952) of the Lebesgue integral.
#'
#' The objective is the negative log-likelihood:
#' \deqn{
#'   -\log\mathcal{L} =
#'   \frac{1}{2}\sum_i q_i
#'   + \sum_i \log\hat{g}(\mathbf{x}_i)
#'   + n\cdot\mathrm{logsumexp}_j\!\left[
#'     -\tfrac{1}{2}q_j - \log\hat{g}(\mathbf{y}_j)
#'   \right]
#' }
#'
#' The \code{"weighted"} model (Jimenez & Soberon 2022, Eq. 5/8) uses
#' \eqn{w = \hat{g}} instead, which is the correct likelihood under the
#' alternative use-availability DGP where presence intensity is
#' proportional to \eqn{f \cdot g}.  See
#' \code{\link{loglik_niche_math_ip_weighted}} for a detailed comparison
#' of both DGP assumptions and known limitations.
#'
#' This function was formerly called \code{niche_kde_bias_corrected}.
#'
#' KDE weights must be precomputed by the user and passed via
#' \code{precomp_w_den}. No KDE is recomputed inside the optimiser.
#'
#' Supports both single-start (\code{start} is a numeric vector) and
#' multi-start (\code{start} is a list of numeric vectors) optimisation.
#'
#' @param occ Numeric matrix (\eqn{n \times p}) of environmental values at
#'   presence points.
#' @param M Numeric matrix (\eqn{n_m \times p}) of environmental values
#'   from the accessible area \eqn{M}. Must have the same columns as
#'   \code{occ}.
#' @param den_idx Integer vector of 1-based row indices selecting the
#'   denominator subset. Must have the same length as \code{precomp_w_den}.
#' @param kde_idx Integer vector of 1-based row indices selecting the KDE
#'   reference subset.
#' @param precomp_w_den Numeric vector of precomputed KDE weights matching
#'   \code{den_idx} in length.
#' @param eta Numeric scalar, shape parameter for the LKJ C-vine prior on
#'   the correlation matrix (default 1 = uniform over correlations).
#' @param start A numeric vector (single-start) or list of numeric vectors
#'   (multi-start). Use [start_theta()] or [start_theta_multiple()] to
#'   generate starting values.
#' @param ... Ignored. Present for wrapper compatibility.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{theta}}{Best parameter vector.}
#'   \item{\code{value}}{Negative log-likelihood at the optimum.}
#'   \item{\code{conv}}{Convergence code (0 = success).}
#'   \item{\code{all_results}}{Data frame of all starts (multi-start only).}
#' }
#'
#' @references
#' Jimenez, L., Soberon, J., Christen, J. A., & Soto, D. (2019).
#' On the problem of modeling a fundamental niche from occurrence data.
#' \emph{Ecological Modelling}, 397, 109823.
#'
#' Jimenez, L., & Soberon, J. (2022).
#' Weighted-normal model for the fundamental niche.
#' \emph{Ecological Modelling}, 438, 109982.
#'
#' Horvitz, D. G., & Thompson, D. J. (1952).
#' A generalization of sampling without replacement from a finite universe.
#' \emph{J. Amer. Statist. Assoc.}, 47(260), 663--685.
#'
#' @seealso [optimize_niche()] for the unified fitting interface,
#'   [loglik_niche_math_ip_weighted()] for the R-level log-likelihood
#'   (includes full DGP derivation and known limitations).
#'
#' @examples
#' \donttest{
#' occ <- as.matrix(example_env_occ_3d)
#' M   <- as.matrix(example_env_m_3d)
#' set.seed(1)
#' den_idx <- sample.int(nrow(M), 300L)
#' kde_idx <- sample.int(nrow(M), 600L)
#' w_den   <- kde_gaussian(M[den_idx, ], M[kde_idx, ])
#' theta0  <- start_theta(example_env_occ_3d)
#' res <- niche_ip_weighted(occ, M, den_idx, kde_idx, w_den, start = theta0)
#' res$value
#' }
#'
#' @export
niche_ip_weighted <- function(occ, M, den_idx, kde_idx, precomp_w_den,
                           eta = 1, start = NULL, ...) {
  # Validate precomputed denominators
  if (length(den_idx) != length(precomp_w_den)) {
    stop("`precomp_w_den` must have the same length as `den_idx` ",
         "(got length(precomp_w_den) = ", length(precomp_w_den),
         ", length(den_idx) = ", length(den_idx), ").")
  }
  if (is.numeric(start) && !all(is.finite(start))) {
    stop("`start` must contain only finite values (no NA, NaN, Inf).")
  }

  # Create compiled C++ objective function
  xptr <- create_niche_obj_ptr(
    env_occ        = as.matrix(occ),
    env_m          = as.matrix(M),
    likelihood     = "ip_weighted",
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
