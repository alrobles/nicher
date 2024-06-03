#' get_suitability
#' Function that calculates the log(suitability)
#' @importFrom mvnfast dmvt
#' @param df A data frame with environmental information
#' @param el_pars A list with ellipse parameter. Contains a vector of means and
#' a covariance matrix. Could came from a MLE after optimization or calculated
#'  from only presence data.
#'
#' @return A vector with suitability
#' @export
#'
#' @examples
#' pars <- get_negloglike_optimr_par(head(spOccPnts, 10), samMPts, fastmode = TRUE, itnmax = 1)
#' get_suitability(spOccPnts, pars)
get_suitability <- function(df, el_pars){
  max.val <- mvnfast::dmvt(X = el_pars$mu,
                           mu = el_pars$mu,
                           sigma = el_pars$S,
                           df = (ncol(df) - 1)
  )
  samMPts <-  mvnfast::dmvt(X = df,
                            mu = el_pars$mu,
                            sigma = el_pars$S,
                            df = (ncol(df) - 1)
  )

  exp(log(samMPts) - log(max.val))

}


