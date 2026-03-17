#' Environmental density from presence points
#' @importFrom mvnfast dmvn
#' @param df A data frame with presence points
#' @param env A raster stack with environmental information
#'
#' @return A data frame with environmental information and the probability density for each point
#' @export
#'
#' @examples
#' get_density(spOccPnts)
get_density <- function(df, env){

  pts <- df
  params <- get_ellip_par(pts)
  dens <- mvnfast::dmvn(X = pts,
                        mu = params$mu,
                        sigma = params$S)
  cbind(pts, dens)
}
