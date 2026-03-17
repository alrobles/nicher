#' Plot method for nicher objects
#'
#' Produces a contour plot of the fitted niche suitability surface for each
#' pair of environmental variables.  For a model with \code{p} variables the
#' output is a \eqn{p \times p}{p x p} panel matrix: diagonal panels show the
#' variable name, lower-triangle panels show a filled contour of the Gaussian
#' density with the presence points overlaid, and upper-triangle panels show
#' the same contour from the transposed perspective.
#'
#' @param x An object of class \code{"nicher"}.
#' @param n_grid Integer; number of grid points along each axis used to
#'   evaluate the suitability surface (default 50).
#' @param env_occ Optional data frame or matrix of presence points to overlay
#'   on the contour plots.  When \code{NULL} (default) no points are drawn.
#' @param env_m Optional data frame or matrix of background points to overlay.
#'   When \code{NULL} (default) no background points are drawn.
#' @param ... Further graphical parameters passed to
#'   \code{\link[graphics]{filled.contour}} or \code{\link[graphics]{contour}}.
#' @return \code{x} invisibly.
#' @export
#' @method plot nicher
#'
#' @examples
#' fit <- nicher(spOccPnts, samMPts, model = "presence_only")
#' plot(fit, env_occ = spOccPnts, env_m = samMPts)
plot.nicher <- function(x, n_grid = 50, env_occ = NULL, env_m = NULL, ...) {
  check_nicher(x)

  mu        <- x$bioscale_params$mu
  S         <- x$bioscale_params$S
  p         <- length(mu)
  env_names <- if (!is.null(x$env_names)) x$env_names else names(mu)
  if (is.null(env_names)) env_names <- paste0("var", seq_len(p))

  col_pairs <- utils::combn(seq_len(p), 2, simplify = FALSE)

  # Gaussian suitability for a pair of variables (marginal 2D projection)
  .suit2d <- function(v1, v2, x_grid, y_grid) {
    mu2  <- mu[c(v1, v2)]
    S2   <- S[c(v1, v2), c(v1, v2)]
    S2_inv <- solve(S2)
    outer(x_grid, y_grid, FUN = function(x, y) {
      mapply(function(xi, yi) {
        d <- c(xi, yi) - mu2
        exp(-0.5 * as.numeric(t(d) %*% S2_inv %*% d))
      }, x, y)
    })
  }

  # Compute index helpers identical to plot.ellipse
  triangle_line_lower <- function(k, n) { 1:k + n * k }
  triangle_line_upper <- function(k, n) { (k + 1):n + n * (k - 1) }

  lower_panel_index <- function(dp) {
    Reduce(c, Map(function(k, n) triangle_line_lower(k, n),
                  1:(dp - 1), dp))
  }
  upper_panel_index <- function(dp) {
    Reduce(c, Map(function(k, n) triangle_line_upper(k, n),
                  1:(dp - 1), dp))
  }

  diagSeq        <- seq(1, p^2, p + 1)
  lowerPanelIdx  <- lower_panel_index(p)
  upperPanelIdx  <- upper_panel_index(p)

  op <- graphics::par(mfrow = c(p, p), mar = c(3, 3, 1, 1))
  on.exit(graphics::par(op))

  heat_cols <- grDevices::hcl.colors(64, "YlOrRd", rev = TRUE)

  for (panel_idx in seq_len(p^2)) {
    if (panel_idx %in% diagSeq) {
      # diagonal: label
      m <- match(panel_idx, diagSeq)
      graphics::plot.new()
      graphics::text(0.5, 0.5, env_names[m], cex = 1.6)

    } else if (panel_idx %in% lowerPanelIdx) {
      m    <- match(panel_idx, lowerPanelIdx)
      pair <- col_pairs[[m]]
      v1   <- pair[1]
      v2   <- pair[2]

      # build grid around [mu - 3sd, mu + 3sd]
      sd1 <- sqrt(S[v1, v1])
      sd2 <- sqrt(S[v2, v2])
      x_grid <- seq(mu[v1] - 3 * sd1, mu[v1] + 3 * sd1, length.out = n_grid)
      y_grid <- seq(mu[v2] - 3 * sd2, mu[v2] + 3 * sd2, length.out = n_grid)
      z      <- .suit2d(v1, v2, x_grid, y_grid)

      graphics::image(x_grid, y_grid, z,
                      col  = heat_cols,
                      xlab = env_names[v1],
                      ylab = env_names[v2])
      graphics::contour(x_grid, y_grid, z,
                        levels     = c(0.1, 0.5, 0.9),
                        add        = TRUE,
                        col        = "white",
                        lwd        = 1.2,
                        drawlabels = TRUE)

      if (!is.null(env_m)) {
        env_m_df <- as.data.frame(env_m)
        graphics::points(env_m_df[[v1]], env_m_df[[v2]],
                         pch = 16, cex = 0.3, col = grDevices::adjustcolor("grey40", 0.5))
      }
      if (!is.null(env_occ)) {
        env_occ_df <- as.data.frame(env_occ)
        graphics::points(env_occ_df[[v1]], env_occ_df[[v2]],
                         pch = 16, cex = 0.8, col = "#1a9850")
      }

    } else if (panel_idx %in% upperPanelIdx) {
      m    <- match(panel_idx, upperPanelIdx)
      pair <- col_pairs[[m]]
      v1   <- pair[1]
      v2   <- pair[2]

      sd1 <- sqrt(S[v1, v1])
      sd2 <- sqrt(S[v2, v2])
      x_grid <- seq(mu[v1] - 3 * sd1, mu[v1] + 3 * sd1, length.out = n_grid)
      y_grid <- seq(mu[v2] - 3 * sd2, mu[v2] + 3 * sd2, length.out = n_grid)
      z      <- .suit2d(v1, v2, x_grid, y_grid)

      graphics::image(y_grid, x_grid, t(z),
                      col  = heat_cols,
                      xlab = env_names[v2],
                      ylab = env_names[v1])
      graphics::contour(y_grid, x_grid, t(z),
                        levels     = c(0.1, 0.5, 0.9),
                        add        = TRUE,
                        col        = "white",
                        lwd        = 1.2,
                        drawlabels = TRUE)

      if (!is.null(env_m)) {
        env_m_df <- as.data.frame(env_m)
        graphics::points(env_m_df[[v2]], env_m_df[[v1]],
                         pch = 16, cex = 0.3, col = grDevices::adjustcolor("grey40", 0.5))
      }
      if (!is.null(env_occ)) {
        env_occ_df <- as.data.frame(env_occ)
        graphics::points(env_occ_df[[v2]], env_occ_df[[v1]],
                         pch = 16, cex = 0.8, col = "#1a9850")
      }
    }
  }

  invisible(x)
}
