#' Negative log-likelihood (inverse-probability-weighted normal, math scale)
#'
#' Computes the negative log-likelihood of the inverse-probability-weighted
#' (IPW) normal niche model.  This model uses \eqn{w = 1/\hat{g}} (the
#' inverse of the background KDE) as importance-sampling weights, applying a
#' Horvitz--Thompson correction so that the denominator approximates the
#' Lebesgue integral \eqn{\int f(\mathbf{x})\,d\mathbf{x}} rather than the
#' \eqn{g}-weighted integral \eqn{\int f\,g\,d\mathbf{x}}:
#'
#' \deqn{
#'   -\log\mathcal{L} =
#'   \frac{1}{2}\sum_i q_i
#'   + \sum_i \log\hat{g}(\mathbf{x}_i)
#'   + n\cdot\mathrm{logsumexp}_j\!\left[
#'     -\tfrac{1}{2}q_j - \log\hat{g}(\mathbf{y}_j)
#'   \right]
#' }
#'
#' where \eqn{q_i} is the Mahalanobis distance at occurrence point
#' \eqn{\mathbf{x}_i}, and the sum over \eqn{j} runs over background
#' points \eqn{\mathbf{y}_j \in M}.
#'
#' @section Data-generating process and motivation:
#'
#' The IPW model is designed for the data-generating process (DGP) of
#' Jimenez et al. (2019), which assumes that presence records are drawn
#' from the fundamental niche density \eqn{f(\mathbf{x};\mu,\Sigma)}
#' restricted to the environments that actually exist in the accessible
#' area \eqn{M}:
#'
#' \deqn{
#'   p(\mathbf{x}_i \mid \text{presence in } M)
#'   = \frac{f(\mathbf{x}_i)}{\int_M f(\mathbf{x})\,d\mathbf{x}}.
#' }
#'
#' The computational challenge is approximating the denominator
#' \eqn{\int_M f\,d\mathbf{x}}.  The background sample
#' \eqn{\{\mathbf{y}_j\}} is drawn from geographic grid cells whose
#' density in environmental space is \eqn{g(\mathbf{x})}.  A naive sum
#' \eqn{\frac{1}{K}\sum_j f(\mathbf{y}_j)} therefore estimates
#' \eqn{\int f\,g\,d\mathbf{x}}, not \eqn{\int f\,d\mathbf{x}}.
#' Dividing each term by \eqn{\hat{g}(\mathbf{y}_j)} corrects the
#' sampling bias via importance sampling:
#'
#' \deqn{
#'   \int_M f(\mathbf{x})\,d\mathbf{x}
#'   = \int_M \frac{f(\mathbf{x})}{g(\mathbf{x})}\,g(\mathbf{x})\,d\mathbf{x}
#'   \approx \frac{1}{K}\sum_{j=1}^{K}
#'     \frac{f(\mathbf{y}_j)}{\hat{g}(\mathbf{y}_j)}.
#' }
#'
#' This is the classical Horvitz--Thompson estimator (Horvitz & Thompson
#' 1952) applied to niche modelling.  The resulting likelihood estimates
#' the fundamental niche on uniform environmental space (Lebesgue measure),
#' removing the distortion that arises when common environments in \eqn{M}
#' dominate the denominator.
#'
#' @section Comparison with the \code{"weighted"} model:
#'
#' The \code{"weighted"} model (Jimenez & Soberon 2022, Eq. 5/8) uses
#' \eqn{w = \hat{g}}, which is the correct likelihood under an alternative
#' DGP where presence intensity is proportional to \eqn{f \cdot g} (the
#' use-availability or resource-selection function model; Warton &
#' Shepherd 2010).  The two models answer different ecological questions:
#'
#' \tabular{lll}{
#'   \strong{Model} \tab \strong{DGP assumption} \tab \strong{Estimates} \cr
#'   \code{"weighted"} \tab \eqn{\lambda \propto f \cdot g} \tab
#'     Niche under M-biased observation process \cr
#'   \code{"ip_weighted"} \tab \eqn{f} restricted to \eqn{M} \tab
#'     Fundamental niche on uniform E-space \cr
#' }
#'
#' Neither model is universally superior.  \code{"weighted"} is appropriate
#' when observation probability correlates with environmental density (e.g.
#' opportunistic citizen-science records).  \code{"ip_weighted"} is
#' appropriate when sampling effort is approximately uniform across \eqn{M}
#' and the goal is to recover the species' intrinsic tolerances independent
#' of which environments happen to be common.
#'
#' @section Known limitations:
#'
#' \enumerate{
#'   \item \strong{Importance-sampling variance.}
#'     When \eqn{\hat{g}(\mathbf{y}_j) \approx 0} (e.g. between
#'     disconnected environmental patches in a multimodal \eqn{M}),
#'     the weights \eqn{1/\hat{g}} can be very large and a handful of
#'     background points may dominate the denominator.  The effective
#'     sample size (ESS) should be monitored; low ESS indicates
#'     unreliable estimates.
#'   \item \strong{KDE bandwidth dependence.}
#'     \eqn{\hat{g}} is computed once using Scott's rule bandwidth.
#'     Over- or under-smoothing can distort the weights.  An adaptive
#'     or cross-validated bandwidth would improve robustness.
#'   \item \strong{Self-regularisation.}
#'     The IPW denominator diverges when \eqn{\sigma \to \infty}
#'     because rare-environment cells contribute \eqn{1/\hat{g} \gg 1}.
#'     This acts as a built-in penalty against unrealistically broad
#'     niches (no explicit ridge prior is needed), but the effective
#'     constraint depends on the tail behaviour of \eqn{\hat{g}} and
#'     can be noisy.
#' }
#'
#' Accepts explicit subsampling indices and precomputed KDE weights for
#' high-performance workflows. Internally delegates to
#' \code{\link{loglik_niche_math_ip_weighted_integrated}}.
#'
#' This function was formerly called
#' \code{loglik_niche_math_kde_bias_corrected}.
#'
#' @param theta Numeric vector of unconstrained parameters (math scale):
#'   \eqn{[\mu, \log\sigma, v]} of length \eqn{2p + p(p-1)/2}.
#' @param env_occ Data frame or matrix (\eqn{n \times p}) of environmental
#'   values at presence points.
#' @param env_m Data frame or matrix (\eqn{n_m \times p}) of background
#'   environmental values from the accessible area M.
#' @param eta Numeric scalar, shape parameter for the LKJ C-vine prior on the
#'   correlation matrix (default 1 = uniform over correlations).
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
#' Warton, D. I., & Shepherd, L. C. (2010).
#' Poisson point process models solve the ``pseudo-absence problem'' for
#' presence-only data in ecology.
#' \emph{Ann. Appl. Stat.}, 4(3), 1383--1402.
#'
#' @seealso [optimize_niche()] for multi-start fitting,
#'   [loglik_niche_math_cpp()] for the M-restricted Gaussian model.
#'
#' @export
#'
#' @examples
#' \donttest{
#' den_idx <- sample.int(nrow(example_env_m_2d), 2000)
#' kde_idx <- sample.int(nrow(example_env_m_2d), 5000)
#' pre_w <- kde_gaussian(
#'   example_env_m_2d[den_idx, ],
#'   example_env_m_2d[kde_idx, ]
#' )
#'
#' loglik_niche_math_ip_weighted(
#'   theta   = start_theta(example_env_occ_2d),
#'   env_occ = example_env_occ_2d,
#'   env_m = example_env_m_2d,
#'   den_idx = den_idx,
#'   kde_idx = kde_idx,
#'   precomp_w_den = pre_w
#' )
#' }
loglik_niche_math_ip_weighted <- function(
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
      loglik_niche_math_ip_weighted_integrated(
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
  loglik_niche_math_ip_weighted_integrated(
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
