#' Predict method for "ellipse" class
#'
#' @param object An object of class "ellipse"
#' @param layers An object of class "rast" from terra package where to evaluate the ellipse model
#' @param returnTable logical if true return table instead a raster with the prediction
#' @param ... Other arguments passed to or from other methods
#' @export predict.ellipse
#' @export
predict.ellipse <- function(object, layers, returnTable = FALSE, ...) {
  rastTab <- terra::as.data.frame(layers, xy = TRUE)
  envData <- rastTab[ ,-c(1:2)]
  rastTabPred <- cbind(rastTab[ ,c(1,2)], suitability = get_suitability(as.matrix(envData), object) )
  rastPred <- terra::rast(rastTabPred)
  if(returnTable == TRUE){
    return(rastTabPred)
  } else {
    return(rastPred)
  }

}
