#' Calculates the inverse of covariance matrix from environmental variables
#' @import stats
#' @param df A matrix or data frame with environmental variables
#'
#' @return A covariance matrix
#' @export
#'
#' @examples
#' get_A_matrix(spOccPnts)
get_A_matrix <-  function(df){
  Sig.ini <- stats::cov(df)
  A.ini <- chol2inv(chol(Sig.ini))
  if(!is.null(colnames(Sig.ini))){
    colnames(A.ini) <- colnames(Sig.ini)
  }
  A.ini
}
