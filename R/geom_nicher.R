# R/geom_nicher.R
# Composable ggplot2 layer constructors for nicher geometry in
# environmental space. Each constructor:
#   * is self-contained: it brings its own data and aes()
#   * references columns positionally via .data[[1]] / .data[[2]]
#   * returns a single ggplot2 layer (composable with `+`)
#   * never reinitializes a ggplot object internally

#' Background environment layer
#'
#' A self-contained \code{ggplot2} layer that draws the environmental
#' background \code{env_m} as a faint point cloud in 2-D environmental
#' space. Designed to be composed with \code{\link{geom_nicher_occ}}
#' and \code{\link{geom_nicher_ellipse}} via the \code{+} operator.
#'
#' @param env_m  Matrix or data frame with at least two columns. Only the
#'   first two columns are used (positional indexing).
#' @param var_names Optional character vector of length 2 used to name the
#'   columns of the layer's internal data frame. Defaults to
#'   \code{c("x1", "x2")}.
#' @param size,alpha,colour Aesthetic parameters; defaults are tuned for a
#'   discreet background layer.
#' @param ... Additional fixed parameters passed to the underlying
#'   \code{ggplot2::geom_point}.
#' @return A \code{ggplot2} layer.
#'
#' @examples
#' \dontrun{
#'   library(ggplot2)
#'   ggplot() + geom_nicher_background(example_env_m_2d)
#' }
#' @export
geom_nicher_background <- function(env_m,
                                    var_names = NULL,
                                    size = 0.3,
                                    alpha = 0.25,
                                    colour = "grey60",
                                    ...) {
  .require_ggplot2("geom_nicher_background")
  d <- .coerce_xy(env_m, var_names, name = "env_m")
  nm <- names(d)
  ggplot2::layer(
    data        = d,
    mapping     = ggplot2::aes(x = .data[[nm[1]]], y = .data[[nm[2]]]),
    geom        = "point",
    stat        = "identity",
    position    = "identity",
    inherit.aes = FALSE,
    show.legend = FALSE,
    params      = list(size = size, alpha = alpha, colour = colour, ...)
  )
}

#' Occurrence-points layer
#'
#' A self-contained \code{ggplot2} layer that draws an occurrence cloud
#' \code{env_occ} in 2-D environmental space.
#'
#' @param env_occ Matrix or data frame with at least two columns. Only the
#'   first two columns are used.
#' @param var_names Optional character vector of length 2 used to name the
#'   columns of the layer's internal data frame.
#' @param size,alpha,colour,shape Aesthetic parameters.
#' @param ... Additional fixed parameters passed to \code{geom_point}.
#' @return A \code{ggplot2} layer.
#'
#' @examples
#' \dontrun{
#'   library(ggplot2)
#'   ggplot() +
#'     geom_nicher_background(example_env_m_2d) +
#'     geom_nicher_occ(example_env_occ_2d)
#' }
#' @export
geom_nicher_occ <- function(env_occ,
                             var_names = NULL,
                             size = 1.2,
                             alpha = 1,
                             colour = "black",
                             shape = 16,
                             ...) {
  .require_ggplot2("geom_nicher_occ")
  d <- .coerce_xy(env_occ, var_names, name = "env_occ")
  nm <- names(d)
  ggplot2::layer(
    data        = d,
    mapping     = ggplot2::aes(x = .data[[nm[1]]], y = .data[[nm[2]]]),
    geom        = "point",
    stat        = "identity",
    position    = "identity",
    inherit.aes = FALSE,
    show.legend = FALSE,
    params      = list(size = size, alpha = alpha,
                       colour = colour, shape = shape, ...)
  )
}

#' Niche-ellipse layer derived from a fitted \code{nicher} model
#'
#' Builds a self-contained \code{ggplot2} layer that traces one or more
#' iso-suitability ellipses implied by a fitted 2-D \code{nicher} model.
#' The ellipse for level \eqn{\alpha} is the contour
#' \eqn{S(x) = \alpha}, equivalently
#' \eqn{(x - \mu)^\top \Sigma^{-1} (x - \mu) = -2 \log \alpha} when
#' \code{level_type = "suitability"}, or the chi-squared confidence ellipse
#' \eqn{(x - \mu)^\top \Sigma^{-1} (x - \mu) =
#'        \mathtt{qchisq}(\alpha, df = 2)} when \code{level_type = "chisq"}.
#'
#' For \code{length(level) > 1} the layer carries one column \code{level}
#' grouping the closed ellipse paths so a single \code{geom_path} draws
#' them all. Map e.g. \code{linetype = factor(level)} downstream to
#' distinguish them visually.
#'
#' @param model A \code{nicher} object with \code{length(model$var_names)
#'   == 2} (or, for legacy fits without \code{var_names}, a 2-D fit).
#' @param level Numeric vector of contour levels in \eqn{(0, 1]}.
#'   Default \code{c(0.95, 0.5, 0.05)} (paper-aligned: core / common /
#'   edge iso-suitability contours).
#' @param level_type Either \code{"suitability"} (the default) or
#'   \code{"chisq"}. See Details.
#' @param n Integer number of points around each ellipse. Default 200.
#' @param linewidth Path linewidth (or \code{size} on ggplot2 < 3.4.0).
#' @param colour Path colour.
#' @param ... Additional fixed parameters passed to \code{geom_path}.
#' @return A \code{ggplot2} layer.
#'
#' @examples
#' \dontrun{
#'   library(ggplot2)
#'   fit <- optimize_niche(
#'     env_occ = example_env_occ_2d,
#'     env_m   = example_env_m_2d,
#'     num_starts = 10L
#'   )
#'   ggplot() +
#'     geom_nicher_background(example_env_m_2d) +
#'     geom_nicher_occ(example_env_occ_2d) +
#'     geom_nicher_ellipse(fit)
#' }
#' @export
geom_nicher_ellipse <- function(model,
                                 level = c(0.95, 0.5, 0.05),
                                 level_type = c("suitability", "chisq"),
                                 n = 200L,
                                 linewidth = 0.6,
                                 colour = "firebrick",
                                 ...) {
  .require_ggplot2("geom_nicher_ellipse")
  .assert_nicher_2d(model, "model")
  level_type <- match.arg(level_type)
  ms <- .recover_mu_sigma(model)
  d  <- .ellipse_paths(ms$mu, ms$Sigma, level, level_type, n,
                       var_names = ms$var_names)
  nm <- names(d)[1:2]

  lw <- .linewidth_param(linewidth)
  ggplot2::layer(
    data        = d,
    mapping     = ggplot2::aes(
                    x = .data[[nm[1]]], y = .data[[nm[2]]],
                    group = .data[["level"]]),
    geom        = "path",
    stat        = "identity",
    position    = "identity",
    inherit.aes = FALSE,
    show.legend = FALSE,
    params      = c(list(colour = colour, ...), lw)
  )
}

#' Iso-suitability contour layer derived from a fitted \code{nicher} model
#'
#' Builds a self-contained \code{ggplot2} layer that draws iso-suitability
#' contours, i.e. level sets \eqn{S(x) = c} of the fitted niche
#' suitability function over a 2-D environmental grid. Unlike
#' \code{\link{geom_nicher_ellipse}}, which traces analytical ellipses
#' and is therefore tied to the Gaussian (multivariate normal) niche
#' geometry, this layer evaluates the model's \emph{actual} suitability
#' function on a grid and contours the result. It works correctly for
#' \strong{any} likelihood family supported by \code{optimize_niche()},
#' including the skew-normal families (\code{"skew_normal"} and
#' \code{"skew_normal_weighted"}) where the iso-suitability sets are
#' \emph{not} ellipses.
#'
#' @section Why "iso-suitability" and not just "ellipse":
#'   For a Gaussian niche, the suitability function
#'   \eqn{S(x) \propto \exp(-\tfrac{1}{2} (x-\mu)^\top \Sigma^{-1}
#'   (x-\mu))} has elliptical level sets; the family of contours
#'   \eqn{S(x) = c} traces nested concentric ellipses around \eqn{\mu}.
#'   For a skew-normal niche
#'   \eqn{S(x) \propto \phi_2(x; \mu, \Sigma) \, \Phi(\alpha^\top
#'   \omega^{-1} (x - \mu))} the \eqn{\Phi(\cdot)} factor breaks the
#'   ellipse symmetry: contours bunch on the side that \eqn{\alpha}
#'   pulls suitability toward, and stretch on the opposite side. The
#'   contour at level \eqn{c} is still a closed curve enclosing
#'   high-suitability environments, but it is no longer an ellipse and
#'   has no closed-form parameterisation. The only honest visualisation
#'   is to evaluate \eqn{S(x)} on a dense 2-D grid and contour the
#'   result -- which is exactly what this layer does.
#'
#' @section Suitability normalisation:
#'   Suitability is always normalised so that \eqn{S(\mu^*) = 1} at the
#'   modal centre: for the Gaussian families that is the centre
#'   \eqn{\mu}; for the skew-normal it is the location parameter
#'   \eqn{\mu} from the SN \eqn{(\mu, \Sigma, \alpha)} parameterisation
#'   (Azzalini & Capitanio 1999), which is generally \emph{not} the
#'   global suitability maximum but is a well-defined reference point.
#'   The contour level \code{level = 0.5} therefore always means
#'   "regions where the niche is at least 50\% as suitable as the
#'   reference centre".
#'
#' @param model A \code{nicher} object with \code{length(model$var_names)
#'   == 2} (or, for legacy fits without \code{var_names}, a 2-D fit).
#' @param level Numeric vector of contour levels in \eqn{(0, 1]}.
#'   Default \code{c(0.95, 0.5, 0.05)}.
#' @param n Integer grid resolution per axis. Default 121 (i.e. a
#'   121 x 121 grid). Higher gives smoother contours at proportionally
#'   higher cost.
#' @param expand Numeric scalar in \eqn{(0, 1)}. The grid spans
#'   \eqn{\mu \pm (1 + \mathtt{expand}) \cdot k \cdot \sigma_{\text{eff}}}
#'   along each axis, where \eqn{k = 4} matches the typical 4-sigma
#'   support of a Gaussian niche and \eqn{\sigma_{\text{eff}}} is the
#'   marginal standard deviation. Default \code{0.1} (10\% padding).
#' @param linewidth Path linewidth (or \code{size} on ggplot2 < 3.4.0).
#' @param colour Path colour.
#' @param ... Additional fixed parameters passed to
#'   \code{\link[ggplot2]{geom_contour}}.
#'
#' @return A \code{ggplot2} layer.
#'
#' @references
#' Azzalini, A. & Capitanio, A. (1999). Statistical applications of the
#' multivariate skew normal distribution.
#' \emph{Journal of the Royal Statistical Society, Series B},
#' \bold{61}(3), 579--602.
#'
#' @examples
#' \dontrun{
#'   library(ggplot2)
#'   fit <- optimize_niche(
#'     env_occ    = example_env_occ_2d,
#'     env_m      = example_env_m_2d,
#'     num_starts = 10L,
#'     likelihood = "skew_normal_weighted"
#'   )
#'   ggplot() +
#'     geom_nicher_background(example_env_m_2d) +
#'     geom_nicher_occ(example_env_occ_2d) +
#'     geom_nicher_isosuitability(fit)
#' }
#' @export
geom_nicher_isosuitability <- function(model,
                                        level = c(0.95, 0.5, 0.05),
                                        n = 121L,
                                        expand = 0.1,
                                        linewidth = 0.6,
                                        colour = "firebrick",
                                        ...) {
  .require_ggplot2("geom_nicher_isosuitability")
  .assert_nicher_2d(model, "model")
  if (!is.numeric(level) || any(!is.finite(level)) ||
      any(level <= 0) || any(level > 1)) {
    stop("`level` must be a numeric vector in (0, 1].", call. = FALSE)
  }
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n < 4L) {
    stop("`n` must be a single integer >= 4.", call. = FALSE)
  }
  if (!is.numeric(expand) || length(expand) != 1L ||
      !is.finite(expand) || expand < 0 || expand >= 1) {
    stop("`expand` must be a single number in [0, 1).", call. = FALSE)
  }

  ms <- .recover_mu_sigma(model)

  # Build a (mu, sigma_eff)-centred grid covering ~ 4 marginal sigmas plus
  # padding. For SN niches the suitability mode is generally offset from
  # mu, but the SN distribution's 99.x% mass still lives within ~4 sigma
  # of mu, so this grid captures the level sets we care about. expand>0
  # adds breathing room so contours don't clip the plotting bbox.
  sigma_eff <- sqrt(diag(ms$Sigma))
  k_sigma   <- 4.0 * (1.0 + expand)
  x1 <- seq(ms$mu[1L] - k_sigma * sigma_eff[1L],
            ms$mu[1L] + k_sigma * sigma_eff[1L], length.out = n)
  x2 <- seq(ms$mu[2L] - k_sigma * sigma_eff[2L],
            ms$mu[2L] + k_sigma * sigma_eff[2L], length.out = n)
  grid_x <- expand.grid(x1 = x1, x2 = x2, KEEP.OUT.ATTRS = FALSE)

  # Evaluate suitability on the grid. We do this in pure R because the
  # grid is small (~10^4 points) and avoiding the C++ blockwise kernel
  # keeps this layer self-contained (it must work in all configurations
  # the model was fitted under).
  diff <- as.matrix(grid_x) - matrix(ms$mu, nrow = nrow(grid_x),
                                     ncol = 2L, byrow = TRUE)
  Sinv <- solve(ms$Sigma)
  q    <- rowSums((diff %*% Sinv) * diff)

  if (ms$is_skew_t) {
    # Skew-t suitability (NCST_p):
    #   S(t) propto (2/A)^{(p+r)/2} * I(t),  A = 1 + q/r
    #   I(t)  = integral_0^inf u^{(p+r)/2 - 1} e^{-u}
    #             * Phi(z(t) sqrt(2 u / (A r))) du
    # Normalised so S(mu) = 1: at t = mu, q = 0, A = 1, z = 0, so
    #   I(mu) = 0.5 * Gamma((p+r)/2),  numerator constants cancel.
    # Approximated by 20-node Gauss-Laguerre (ample for 2D grid plots;
    # the expensive C++ kernel uses 32 nodes for 6-D fits).
    r   <- exp(ms$log_r)
    pdim <- 2L
    a   <- (pdim + r) / 2.0
    z   <- as.numeric(diff %*% (ms$alpha / sigma_eff))
    A   <- 1.0 + q / r
    # Hard-coded 20-node standard Gauss-Laguerre nodes/weights (laggauss(20)).
    gl <- .gl_quad_20()
    # For each grid row evaluate the integral via vectorised matrix product.
    # Each column corresponds to a quadrature node.
    log_u <- log(gl$nodes)
    # u^{a-1} -> exp((a-1) log u). Build matrix [n_grid x Q] of arguments.
    pnorm_arg <- outer(z / sqrt(A), sqrt(2 * gl$nodes / r))
    # Integrand: u^{a-1} * Phi(z sqrt(2u/(A r))). The exp(-u) is folded
    # into the Gauss-Laguerre weights already.
    integrand <- exp((a - 1.0) * log_u)[col(pnorm_arg)] *
                   stats::pnorm(pnorm_arg)
    integral  <- as.numeric(integrand %*% gl$weights)
    # I(mu) = 0.5 * Gamma(a) -> normalising constant.
    I_mu <- 0.5 * gamma(a)
    S    <- (2 / A)^a * integral / ((2)^a * I_mu)
  } else if (ms$is_skew) {
    # Skew-normal suitability: phi_2(.; mu, Sigma) * Phi(alpha . omega^-1 d)
    # divided by its value at mu (where the Gaussian factor = 1 and
    # Phi(0) = 0.5), so S(mu) = 1 by construction.
    z   <- as.numeric(diff %*% (ms$alpha / sigma_eff))
    S   <- exp(-0.5 * q) * stats::pnorm(z) / 0.5
  } else {
    S <- exp(-0.5 * q)
  }
  # Numerical safety: clip [0, 1] (small overshoot can occur when SN
  # mode > mu, but contours of interest live within (0, 1)).
  S[!is.finite(S)] <- 0
  S <- pmin(pmax(S, 0), 1)

  d <- data.frame(grid_x, S = S)
  vn <- ms$var_names
  if (!is.null(vn) && length(vn) == 2L) {
    names(d)[1:2] <- vn
  }
  nm <- names(d)[1:2]

  lw <- .linewidth_param(linewidth)
  ggplot2::layer(
    data        = d,
    mapping     = ggplot2::aes(
                    x = .data[[nm[1]]], y = .data[[nm[2]]],
                    z = .data[["S"]]),
    geom        = "contour",
    stat        = "contour",
    position    = "identity",
    inherit.aes = FALSE,
    show.legend = FALSE,
    params      = c(list(breaks = sort(unique(level), decreasing = TRUE),
                         colour = colour, ...), lw)
  )
}
