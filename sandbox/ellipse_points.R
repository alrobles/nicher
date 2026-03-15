#' Ellipse points
#'
#' @param par_ellipse An ellipse object with ellipse parameters
#' @param columns The number of columns of the ellipse to get the points
#' @param n_points The number of points for the ellipses to plot
#'
#' @return A data frame with the coordinates to plot
#' @export
#'
#' @examples
#' ellipse_points(ellipse_par_example)
#' ellipse_points(ellipse_par_example, columns = c(1, 3)  )
ellipse_points <- function(par_ellipse, columns = c(1, 2), n_points = 100){

  theta <- par_ellipse[[1]][columns]
  sigma <- par_ellipse[[2]][columns, columns]

  # we start from points on the unit circle
  n_points

  xy <- cbind(sin(seq(0, 2 * pi, length.out = n_points)),
              cos(seq(0, 2 * pi, length.out = n_points)))

  # then we scale the dimensions
  ev <- eigen(sigma)
  xy[, 1] <- xy[, 1] * 1
  xy[, 2] <- xy[, 2] * sqrt(min(ev$values) / max(ev$values))

  # then rotate
  phi <- atan(ev$vectors[2, 1] / ev$vectors[1, 1])
  R <- matrix(c(cos(phi), sin(phi), -sin(phi), cos(phi)), 2)
  xy <- tcrossprod(R, xy)

  # the quantiles you ask for
  chi_vals <- qchisq(c(0.95), df = 2) * max(ev$values)

  #####
  # Plot contours
  df <- data.frame(x = sqrt(chi_vals) * xy[1, ] + theta[1] ,
                   y  = sqrt(chi_vals) * xy[2, ] + theta[2])
  df$variable <-  paste0(columns, collapse = "_")
  return(df)
}
