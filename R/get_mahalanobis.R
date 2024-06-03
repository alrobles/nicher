#' Mahalanobis distance
#'
#' @param df A data frame with environmental information according to presence
#'   points
#' @param el_pars A list with ellipsoid parameters
#'
#' @return A vector with Mahalanobis distance according to each point in
#'   environmental space
#' @export
#'
#' @examples
#' pars <- get_negloglike_optimr_par(spOccPnts, samMPts, lower = TRUE, fastmode = TRUE )
#' get_mahalanobis(spOccPnts, pars)
get_mahalanobis <- function(df, el_pars){
  mahadist <- mvnfast::maha(X = df,
                            mu = el_pars$mu,
                            sigma = el_pars$S)
  return(mahadist)
}
