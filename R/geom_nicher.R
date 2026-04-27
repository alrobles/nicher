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
